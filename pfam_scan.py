#!/usr/bin/env python3
"""
pfam_qc.py – build a Pfam presence/absence matrix from a directory
              of HMMER *.domtblout (or *.out) files.

Changes relative to the original:
  * Improved overlap resolver (protein-coords, 15 aa + 25 % rule).
  * Bitscore-based priority.
  * Optional parallel processing with --cpu and tqdm progress bar.
  * User-tunable overlap thresholds (--ov-min-len / --ov-min-frac).
  * Revised verbosity: -v for INFO, -vv for DEBUG, default is WARNING (suppress INFO).
  * Fixed critical bugs and improved error handling.
  * Added support for both hmmscan and hmmsearch formats with --mode flag.
  * ADDED: GA threshold support for Pfam-style family-specific thresholds.
"""
from __future__ import annotations
import os, glob, sys, argparse, logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from typing import List, Dict, Tuple, Optional
import numpy as np
import pandas as pd
from tqdm import tqdm

# ---------------------------------------------------------------------
# Column mappings for HMMER domtblout format
# ---------------------------------------------------------------------
# HMMSCAN: searches sequences against HMM database (e.g., proteins vs Pfam)
# target = HMM (Pfam), query = sequence (protein)
HMMSCAN_COLS = {
    "pfam_full" : 1,     # target accession (PFxxxxx.version)
    "protein"   : 3,     # query name (protein ID)
    "hmm_len"   : 2,     # target length (HMM length)
    "bitscore"  : 7,     # full sequence score
    "iEvalue"   : 12,    # independent E-value
    "hmm_from"  : 15,    # hmm coord from
    "hmm_to"    : 16,    # hmm coord to
    "ali_from"  : 17,    # ali coord from
    "ali_to"    : 18,    # ali coord to
}

# HMMSEARCH: searches HMM against sequence database (e.g., Pfam vs proteins)
# target = sequence (protein), query = HMM (Pfam)
HMMSEARCH_COLS = {
    "protein"   : 0,     # target name (protein ID)
    "pfam_full" : 3,     # query name (Pfam ID)
    "hmm_len"   : 5,     # query length (HMM length)
    "bitscore"  : 7,     # full sequence score
    "iEvalue"   : 6,     # full sequence E-value
    "hmm_from"  : 15,    # hmm coord from
    "hmm_to"    : 16,    # hmm coord to
    "ali_from"  : 17,    # ali coord from
    "ali_to"    : 18,    # ali coord to
}

# ---------------------------------------------------------------------
# GA threshold parser
# ---------------------------------------------------------------------
def parse_ga_thresholds(hmm_file: str) -> Dict[str, Dict[str, float]]:
    """
    Parse GA (Gathering) thresholds from Pfam-A.hmm file.
    
    Returns:
        Dict mapping Pfam ID to {'domain_ga': float, 'seq_ga': float}
    """
    ga_thresholds = {}
    
    if not os.path.exists(hmm_file):
        logging.warning(f"GA threshold file not found: {hmm_file}")
        return ga_thresholds
    
    logging.info(f"Parsing GA thresholds from {hmm_file}")
    
    try:
        with open(hmm_file, 'r', encoding='utf-8') as f:
            current_acc = None
            line_count = 0
            
            for line in f:
                line_count += 1
                
                # Progress indicator for large files
                if line_count % 100000 == 0:
                    logging.debug(f"Processed {line_count:,} lines, found {len(ga_thresholds):,} GA thresholds")
                
                if line.startswith('ACC   '):
                    # Extract Pfam accession (remove version if present)
                    try:
                        acc_line = line.strip().split(None, 1)[1]  # Get everything after ACC
                        current_acc = acc_line.split('.')[0]  # Remove version (e.g., PF00001.23 -> PF00001)
                    except (IndexError, AttributeError):
                        logging.warning(f"Malformed ACC line at line {line_count}: {line.strip()}")
                        current_acc = None
                        
                elif line.startswith('GA   ') and current_acc:
                    # Parse GA line format: "GA   25.00 25.00;"
                    try:
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            domain_ga = float(parts[1])
                            seq_ga = float(parts[2].rstrip(';'))
                            ga_thresholds[current_acc] = {
                                'domain_ga': domain_ga,
                                'seq_ga': seq_ga
                            }
                        else:
                            logging.warning(f"Malformed GA line at line {line_count}: {line.strip()}")
                    except (ValueError, IndexError) as e:
                        logging.warning(f"Error parsing GA line at line {line_count}: {e}")
                        
                elif line.startswith('//'):
                    # End of record
                    current_acc = None
                    
    except IOError as e:
        logging.error(f"Cannot read GA threshold file {hmm_file}: {e}")
        return {}
    
    logging.info(f"Parsed {len(ga_thresholds):,} GA thresholds")
    return ga_thresholds

def parse_domtbl(domtbl: str, num_cols: int, mode: str) -> List[Dict]:
    """
    Parse one domtblout file -> list of hit-dicts.
    Skips comment / blank lines and lines with < num_cols columns.
    
    Args:
        domtbl: Path to domtblout file
        num_cols: Expected number of columns
        mode: Either 'hmmscan' or 'hmmsearch'
    """
    hits = []
    genome = os.path.basename(domtbl).rsplit(".", 1)[0]
    logging.debug(f"Parsing {domtbl} (⇒ genome '{genome}', mode='{mode}')")
    
    # Select appropriate column mapping
    if mode == "hmmscan":
        cols_map = HMMSCAN_COLS
    elif mode == "hmmsearch":
        cols_map = HMMSEARCH_COLS
    else:
        raise ValueError(f"Invalid mode: {mode}. Must be 'hmmscan' or 'hmmsearch'")
    
    if not os.path.exists(domtbl):
        logging.error(f"File not found: {domtbl}")
        return hits
    
    try:
        with open(domtbl, 'r', encoding='utf-8') as fh:
            for lineno, line in enumerate(fh, 1):
                if line.startswith("#") or not line.strip():
                    continue
                
                cols = line.rstrip().split()
                if len(cols) < num_cols:
                    if lineno <= 5 or lineno % 1000 == 0:  # Show first few errors
                        logging.warning(
                            f"Short line {lineno} in {domtbl!r}: "
                            f"{len(cols)} < {num_cols} columns – skipped"
                        )
                    continue
                
                try:
                    # Extract Pfam ID (remove version if present)
                    pfam_full = cols[cols_map["pfam_full"]]
                    pfam_id = pfam_full.split(".")[0] if "." in pfam_full else pfam_full
                    
                    rec = {
                        "genome"   : genome,
                        "protein"  : cols[cols_map["protein"]],
                        "pfam"     : pfam_id,
                        "bitscore" : float(cols[cols_map["bitscore"]]),
                        "iEvalue"  : float(cols[cols_map["iEvalue"]]),
                        "hmm_from" : int(cols[cols_map["hmm_from"]]),
                        "hmm_to"   : int(cols[cols_map["hmm_to"]]),
                        "ali_from" : int(cols[cols_map["ali_from"]]),
                        "ali_to"   : int(cols[cols_map["ali_to"]]),
                        "hmm_len"  : int(cols[cols_map["hmm_len"]]),
                    }
                    
                    # Validate coordinates
                    if (rec["hmm_from"] > rec["hmm_to"] or 
                        rec["ali_from"] > rec["ali_to"] or
                        rec["hmm_from"] < 1 or rec["ali_from"] < 1):
                        logging.warning(f"Invalid coordinates at line {lineno} in {domtbl}")
                        continue
                        
                    hits.append(rec)
                    
                except (IndexError, ValueError) as e:
                    logging.error(f"Parse error {domtbl}:{lineno}: {e}")
                    continue
                    
    except IOError as e:
        logging.error(f"Cannot read file {domtbl}: {e}")
        return []
    
    logging.info(f"Parsed {len(hits):,} raw hits from {domtbl}")
    return hits

# ---------------------------------------------------------------------
# QC filter with GA threshold support
# ---------------------------------------------------------------------
def qc_filter(
    hits: List[Dict], 
    evalue_thresh: float, 
    cov_thresh: float,
    ga_thresholds: Optional[Dict[str, Dict[str, float]]] = None,
    use_ga: bool = False
) -> List[Dict]:
    """
    Apply E-value + HMM coverage filter, optionally using GA thresholds.
    
    Args:
        hits: List of hit dictionaries
        evalue_thresh: E-value threshold (used as fallback when GA not available)
        cov_thresh: Coverage threshold
        ga_thresholds: Dict of GA thresholds per family
        use_ga: Whether to use GA thresholds when available
    """
    if not hits:
        return []
    
    passed = []
    ga_used = 0
    fallback_used = 0
    
    for h in hits:
        try:
            # Calculate HMM coverage
            hmm_span = h["hmm_to"] - h["hmm_from"] + 1
            cov = hmm_span / h["hmm_len"] if h["hmm_len"] > 0 else 0
            
            # Coverage filter (always applied)
            if cov < cov_thresh:
                continue
                
            # Score/E-value filter
            pass_score = False
            
            if use_ga and ga_thresholds and h["pfam"] in ga_thresholds:
                # Use GA threshold (domain threshold for individual hits)
                domain_ga = ga_thresholds[h["pfam"]]["domain_ga"]
                if h["bitscore"] >= domain_ga:
                    pass_score = True
                    ga_used += 1
                    logging.debug(f"GA filter: {h['pfam']} score={h['bitscore']:.1f} >= {domain_ga:.1f}")
            else:
                # Fallback to E-value threshold
                if h["iEvalue"] <= evalue_thresh:
                    pass_score = True
                    fallback_used += 1
                    if use_ga:
                        logging.debug(f"Fallback filter: {h['pfam']} (no GA available)")
            
            if pass_score:
                h["coverage"] = cov
                passed.append(h)
                
        except (KeyError, ZeroDivisionError) as e:
            logging.warning(f"QC filter error for hit: {e}")
            continue
    
    if use_ga:
        logging.info(f"QC filter: kept {len(passed):,}/{len(hits):,} hits "
                    f"({ga_used:,} via GA, {fallback_used:,} via E-value)")
    else:
        logging.info(f"QC filter: kept {len(passed):,}/{len(hits):,} hits")
    
    return passed

# ---------------------------------------------------------------------
# Improved overlap resolver
# ---------------------------------------------------------------------
def resolve_overlaps(
    hits: List[Dict],
    min_len: int = 15,
    min_frac: float = 0.25
) -> List[Dict]:
    """
    Remove *substantially* overlapping domains on the same protein.

    Conflict if:
      overlap_len > min_len  AND
      overlap_len / shorter_interval ≥ min_frac  (default 15 aa + 25 %).

    Priority: higher bitscore first (ties: lower E-value, then longer alignment).
    """
    if not hits:
        return []
    
    by_prot: Dict[Tuple[str, str], List[Dict]] = {}
    retained = []

    # Bucket by (genome, protein)
    for h in hits:
        key = (h["genome"], h["protein"])
        by_prot.setdefault(key, []).append(h)

    for group in by_prot.values():
        if not group:
            continue
            
        # Sort by priority: descending bitscore, ascending E-value, descending alignment length
        group.sort(
            key=lambda x: (-x["bitscore"], x["iEvalue"], -(x["ali_to"] - x["ali_from"]))
        )
        
        accepted = []
        for h in group:
            clash = False
            h_len = h["ali_to"] - h["ali_from"] + 1
            
            for g in accepted:
                # Calculate overlap in protein coordinates
                overlap_start = max(h["ali_from"], g["ali_from"])
                overlap_end = min(h["ali_to"], g["ali_to"])
                overlap_len = max(0, overlap_end - overlap_start + 1)
                
                if overlap_len <= min_len:
                    continue
                
                # Check if overlap is substantial
                g_len = g["ali_to"] - g["ali_from"] + 1
                shorter_len = min(h_len, g_len)
                
                if shorter_len > 0 and overlap_len / shorter_len >= min_frac:
                    clash = True
                    logging.debug(
                        f"Overlap clash: {h['pfam']} vs {g['pfam']} on "
                        f"{h['genome']}:{h['protein']} ({overlap_len}/{shorter_len})"
                    )
                    break
            
            if not clash:
                accepted.append(h)
        
        retained.extend(accepted)

    logging.info(f"Overlap resolver: {len(retained):,}/{len(hits):,} retained")
    return retained

# ---------------------------------------------------------------------
# Helper for per-file processing (so it can be parallelised)
# ---------------------------------------------------------------------
def process_file(
    path: str,
    num_cols: int,
    evalue: float,
    cov: float,
    do_overlap: bool,
    ov_min_len: int,
    ov_min_frac: float,
    mode: str,
    ga_thresholds: Optional[Dict[str, Dict[str, float]]] = None,
    use_ga: bool = False
) -> List[Dict]:
    """Process a single domtblout file through the full pipeline."""
    try:
        hits = parse_domtbl(path, num_cols, mode)
        if not hits:
            return []
        
        hits = qc_filter(hits, evalue, cov, ga_thresholds, use_ga)
        if not hits:
            return []
        
        if do_overlap:
            hits = resolve_overlaps(hits, ov_min_len, ov_min_frac)
        
        return hits
    except Exception as e:
        logging.error(f"Error processing file {path}: {e}")
        return []

# ---------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(
        description="QC-filtered Pfam presence/absence matrix from domtblout files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--domtbls", required=True,
                   help="Directory of domtblout (or .out) files")
    p.add_argument("--ext", default=".domtblout",
                   help="File extension to scan for")
    p.add_argument("--out", default="pfam_presence.tsv",
                   help="Output TSV file")
    p.add_argument("--mode", choices=["hmmscan", "hmmsearch"], default="hmmscan",
                   help="HMMER mode: hmmscan (seqs vs HMMs) or hmmsearch (HMMs vs seqs)")
    p.add_argument("--num_cols", type=int, default=23,
                   help="Expected columns in domtblout")
    p.add_argument("--evalue", type=float, default=1e-3,
                   help="Max independent E-value threshold (fallback when GA not available)")
    p.add_argument("--cov", type=float, default=0.5,
                   help="Min HMM coverage fraction")
    p.add_argument("--resolve-overlaps", action="store_true",
                   help="Apply Pfam-style overlap resolver")
    p.add_argument("--ov-min-len", type=int, default=15,
                   help="Absolute residues for clash")
    p.add_argument("--ov-min-frac", type=float, default=0.25,
                   help="Fraction of shorter interval")
    p.add_argument("--cpu", type=int, default=1,
                   help="Number of worker processes")
    p.add_argument("-v", "--verbose", action="count", default=0,
                   help="Increase verbosity: -v = INFO, -vv = DEBUG")
    p.add_argument("--per-sequence", action="store_true",
                   help="Build matrix per protein sequence (columns = sequences) "
                        "instead of per file (columns = genomes/files)")
    p.add_argument("--counts", action="store_true",
                   help="Output raw Pfam gene counts instead of binary presence")
    p.add_argument("--log-transform", action="store_true",
                   help="Apply log(count + pseudocount) – requires --counts")
    p.add_argument("--pseudocount", type=float, default=0.5,
                   help="Pseudocount added before log-transform (default 0.5)")
    
    # NEW: GA threshold support
    p.add_argument("--ga-file", type=str, default=None,
                   help="Path to Pfam-A.hmm file for GA threshold parsing")
    p.add_argument("--use-ga", action="store_true",
                   help="Use GA thresholds when available (requires --ga-file)")
    
    args = p.parse_args()

    if args.log_transform and not args.counts:
        p.error("--log-transform requires --counts")
    if args.pseudocount <= 0:
        p.error("--pseudocount must be > 0")
    if args.use_ga and not args.ga_file:
        p.error("--use-ga requires --ga-file")
    
    # Validate arguments
    if args.evalue <= 0:
        logging.error("E-value threshold must be positive")
        sys.exit(1)
    if not 0 <= args.cov <= 1:
        logging.error("Coverage threshold must be between 0 and 1")
        sys.exit(1)
    if args.ov_min_frac < 0 or args.ov_min_frac > 1:
        logging.error("Overlap fraction must be between 0 and 1")
        sys.exit(1)
    if args.cpu < 1:
        logging.error("CPU count must be at least 1")
        sys.exit(1)

    # -----------------------------------------------------------------
    # Set logging levels based on -v count:
    #   0 → WARNING, 1 → INFO, 2+ → DEBUG
    if args.verbose >= 2:
        level = logging.DEBUG
    elif args.verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING
    
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    
    filter_mode = "GA thresholds" if args.use_ga else "E-value"
    logging.info(f"Starting Pfam QC pipeline (mode: {args.mode}, filter: {filter_mode})")

    # -----------------------------------------------------------------
    # Parse GA thresholds if requested
    ga_thresholds = None
    if args.use_ga:
        ga_thresholds = parse_ga_thresholds(args.ga_file)
        if not ga_thresholds:
            logging.error("No GA thresholds found - check your --ga-file path")
            sys.exit(1)

    # -----------------------------------------------------------------
    # Locate files
    if not os.path.isdir(args.domtbls):
        logging.error(f"Directory not found: {args.domtbls}")
        sys.exit(1)
    
    pattern = os.path.join(args.domtbls, f"*{args.ext}")
    files = sorted(glob.glob(pattern))
    
    if not files:
        logging.error(f"No files found matching {pattern}")
        sys.exit(1)
    
    logging.info(f"Found {len(files):,} files with extension '{args.ext}'")

    # -----------------------------------------------------------------
    # Prepare partial function for pool workers
    work = partial(
        process_file,
        num_cols=args.num_cols,
        evalue=args.evalue,
        cov=args.cov,
        do_overlap=args.resolve_overlaps,
        ov_min_len=args.ov_min_len,
        ov_min_frac=args.ov_min_frac,
        mode=args.mode,
        ga_thresholds=ga_thresholds,
        use_ga=args.use_ga
    )

    all_hits: List[Dict] = []

    # -----------------------------------------------------------------
    # Process files
    if args.cpu > 1:
        logging.info(f"Processing in parallel with {args.cpu} CPU workers")
        with ProcessPoolExecutor(max_workers=args.cpu) as pool:
            futures = {pool.submit(work, f): f for f in files}
            try:
                for fut in tqdm(as_completed(futures), total=len(files), unit="file"):
                    result = fut.result()
                    all_hits.extend(result)
            except KeyboardInterrupt:
                logging.warning("Interrupted by user")
                pool.shutdown(wait=False)
                sys.exit(1)
    else:
        logging.info("Processing sequentially (CPU=1)")
        try:
            for f in tqdm(files, unit="file"):
                result = work(f)
                all_hits.extend(result)
        except KeyboardInterrupt:
            logging.warning("Interrupted by user")
            sys.exit(1)

    # -----------------------------------------------------------------
    # Build and save matrix
    if not all_hits:
        logging.error("No QC-passing hits found – nothing to pivot.")
        sys.exit(1)

    logging.info(f"Total hits after QC: {len(all_hits):,}")
    
    try:
        df = pd.DataFrame(all_hits)

        # Column axis: per-file (genome) by default, per-sequence with --per-sequence
        col_axis = "protein" if args.per_sequence else "genome"
        col_label = "sequences" if args.per_sequence else "genomes"

        if args.counts:            # ───────────── raw counts branch ─────────────
            df["count"] = 1
            mat = df.pivot_table(index="pfam",
                                columns=col_axis,
                                values="count",
                                aggfunc="sum",
                                fill_value=0)

            if args.log_transform:                 # log(count + pseudocount)
                mat = np.log(mat + args.pseudocount)

        else:                     # ───────────── binary presence/absence ─────────────
            df["presence"] = 1
            mat = df.pivot_table(index="pfam",
                                columns=col_axis,
                                values="presence",
                                aggfunc="max",
                                fill_value=0)


        # Ensure output directory exists
        out_dir = os.path.dirname(args.out)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        logging.info(f"Matrix dimensions: {mat.shape[0]} Pfam families × {mat.shape[1]} {col_label}")
        
        # Save matrix
        mat.to_csv(args.out, sep="\t")
        
        logging.info(f"Wrote QC-filtered matrix to {args.out}")
        
    except Exception as e:
        logging.error(f"Error building matrix: {e}")
        sys.exit(1)


# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()