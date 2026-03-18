#!/usr/bin/env python3
"""
pfam_scan_mp.py – multiprocessing *and* multithreading Pfam scans
"""

from __future__ import annotations
import argparse, math, os, sys, multiprocessing as mp
from pathlib import Path
from typing import Union, List, Tuple

from tqdm import tqdm
from pyhmmer import easel, plan7, hmmer          # pyhmmer ≥ 0.11.1

# Global variable for worker processes
PROFILES = None

# ───────────────────────── CLI ─────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate *.faa with Pfam-A using N processes × M threads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--pfam",       required=True, type=Path)
    p.add_argument("--faa-dir",    default=Path("."), type=Path)
    p.add_argument("--output-dir", default=None, type=Path)
    p.add_argument("--threads",    default=0, type=int,
                   help="Threads *per* process (0 → all cores / Nproc)")
    p.add_argument("--processes",  default=1, type=int,
                   help="How many processes scan files in parallel")
    p.add_argument("--stream",     action="store_true",
                   help="Stream Pfam from disk (low RAM, slower)")
    return p.parse_args()

# ─────────────────── Pfam loader ───────────────────
def load_pfam(pfam: Path, alphabet: easel.Alphabet,
              stream: bool) -> Union[plan7.OptimizedProfileBlock,
                                     plan7.HMMPressedFile]:
    """Load Pfam database either into memory or as a streamed file."""
    try:
        if stream:
            pressed = plan7.HMMPressedFile(str(pfam))
            print(f"▶ Pfam will be streamed ({len(pressed):,} models).")
            return pressed
        print("▶ Loading Pfam into RAM …")
        block = plan7.OptimizedProfileBlock(alphabet)
        with plan7.HMMFile(str(pfam)) as db:
            for op in tqdm(db.optimized_profiles(), desc="Prefetch",
                           unit="model", colour="cyan"):
                block.append(op)
        print(f"  Cached {len(block):,} models (≈{len(block)*3500/1e6:.1f} MB).")
        return block
    except Exception as e:
        sys.exit(f"❌ Error loading Pfam database: {e}")

# ───────────────── Global initializer for worker processes ─────────────────
def init_worker(profiles_data):
    """Initialize worker process with pre-loaded Pfam database."""
    global PROFILES
    PROFILES = profiles_data

# ───────────────── worker function ─────────────────
def worker_task(args: Tuple[str, str, int, int]) -> Tuple[str, int]:
    """
    Run hmmscan on one file using the globally initialized PROFILES.
    """
    faa, out, threads, Z = args
    
    global PROFILES
    if PROFILES is None:
        raise RuntimeError("Worker not properly initialized - PROFILES is None")
    
    try:
        alphabet = easel.Alphabet.amino()
        
        # Create progress bar for this file
        inner = tqdm(desc=Path(faa).name, unit="seq",
                     colour="green", leave=False, position=0)
        
        def callback(query, hit_num, bar=inner):
            bar.update()
        
        # Process sequences
        with easel.SequenceFile(faa, digital=True, alphabet=alphabet) as sf, \
             open(out, "wb") as fh:  # Keep as binary mode for pyhmmer output
            
            hit_count = 0
            for hits in hmmer.hmmscan(sf, PROFILES,
                                      cpus=threads, Z=Z, callback=callback):
                # Write hits in binary mode as pyhmmer expects
                hits.write(fh, format="domains")
                hit_count += len(hits)
        
        inner.close()
        
        # Return file info and hit count
        file_size = os.path.getsize(out) if os.path.exists(out) else 0
        return (faa, file_size)
        
    except Exception as e:
        if 'inner' in locals():
            inner.close()
        raise RuntimeError(f"Error processing {faa}: {e}")

# ───────────────────────── main ─────────────────────────
def main():
    a = parse_args()
    
    # Validate inputs
    if not a.pfam.exists(): 
        sys.exit(f"❌ Pfam not found: {a.pfam}")
    
    if not a.faa_dir.exists():
        sys.exit(f"❌ Input directory not found: {a.faa_dir}")
    
    # Create output directory if specified
    if a.output_dir: 
        a.output_dir.mkdir(parents=True, exist_ok=True)

    # Set multiprocessing start method
    # Fork is more efficient but not available on all platforms
    try: 
        mp.set_start_method("fork")
        print("▶ Using 'fork' multiprocessing method")
    except RuntimeError: 
        print("▶ Using default multiprocessing method (spawn/forkserver)")

    # Find input files
    faa_files = sorted(a.faa_dir.glob("*.faa"))
    if not faa_files: 
        sys.exit(f"❌ No *.faa files found in {a.faa_dir}")
    
    print(f"▶ Found {len(faa_files)} .faa files to process")

    # Calculate optimal thread distribution
    hw = a.threads or os.cpu_count() or 1
    nproc = max(1, min(a.processes, len(faa_files)))  # Don't use more processes than files
    thr   = max(1, math.ceil(hw / nproc))
    print(f"▶ {nproc} process(es) × {thr} thread(s) each (≤ {hw} total)")

    # Load Pfam in main process for sharing via fork
    try:
        alphabet = easel.Alphabet.amino()
        print("▶ Loading Pfam database for sharing across processes...")
        profiles_main = load_pfam(a.pfam, alphabet, a.stream)
        Z = len(profiles_main)
        print(f"▶ Database contains {Z:,} models")
    except Exception as e:
        sys.exit(f"❌ Failed to load Pfam database: {e}")

    # Prepare tasks
    tasks: List[Tuple[str, str, int, int]] = []
    for faa in faa_files:
        if a.output_dir:
            out = (a.output_dir / faa.name).with_suffix(".domtblout")
        else:
            out = faa.with_suffix(".domtblout")
        tasks.append((str(faa), str(out), thr, Z))

    # Process files
    print(f"▶ Starting annotation of {len(tasks)} files...")
    
    try:
        bar = tqdm(total=len(tasks), desc="Annotating files", unit="file")
        
        # Pass the loaded profiles directly to workers via fork
        with mp.Pool(processes=nproc, 
                     initializer=init_worker,
                     initargs=(profiles_main,)) as pool:
            
            completed = 0
            for result in pool.imap_unordered(worker_task, tasks):
                faa_file, file_size = result
                bar.set_postfix({"last": Path(faa_file).name, 
                                "size": f"{file_size/1024:.1f}KB"})
                bar.update()
                completed += 1
        
        bar.close()
        print(f"✓ Successfully annotated {completed}/{len(tasks)} files.")
        
    except KeyboardInterrupt:
        print("\n❌ Interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error during processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()