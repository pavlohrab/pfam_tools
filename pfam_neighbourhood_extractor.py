#!/usr/bin/env python3
"""
pfam_neighbourhood_extractor.py - Extract genomic neighbourhoods around Pfam domains
from GenBank files based on HMMER domtblout results.

Uses .domtblout files from run_hmmscan.py together with GenBank files.

This script:
1. Parses HMMER domtblout files and applies QC filtering
2. Optionally resolves overlapping domains
3. Reads corresponding GenBank files
4. Extracts genomic neighbourhoods (N genes upstream/downstream) around target Pfam domains
5. Outputs new GenBank files with the neighbourhoods
6. Optionally creates .faa files with protein sequences per Pfam domain
"""

import os
import glob
import sys
import argparse
import logging
import re
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import multiprocessing

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


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
        with open(hmm_file, "r", encoding="utf-8") as f:
            current_acc = None
            line_count = 0

            for line in f:
                line_count += 1

                if line_count % 100000 == 0:
                    logging.debug(
                        f"Processed {line_count:,} lines, "
                        f"found {len(ga_thresholds):,} GA thresholds"
                    )

                if line.startswith("ACC   "):
                    try:
                        acc_line = line.strip().split(None, 1)[1]
                        current_acc = acc_line.split(".")[0]
                    except (IndexError, AttributeError):
                        logging.warning(
                            f"Malformed ACC line at line {line_count}: {line.strip()}"
                        )
                        current_acc = None

                elif line.startswith("GA   ") and current_acc:
                    try:
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            domain_ga = float(parts[1])
                            seq_ga = float(parts[2].rstrip(";"))
                            ga_thresholds[current_acc] = {
                                "domain_ga": domain_ga,
                                "seq_ga": seq_ga,
                            }
                        else:
                            logging.warning(
                                f"Malformed GA line at line {line_count}: {line.strip()}"
                            )
                    except (ValueError, IndexError) as e:
                        logging.warning(
                            f"Error parsing GA line at line {line_count}: {e}"
                        )

                elif line.startswith("//"):
                    current_acc = None

    except IOError as e:
        logging.error(f"Cannot read GA threshold file {hmm_file}: {e}")
        return {}

    logging.info(f"Parsed {len(ga_thresholds):,} GA thresholds")
    return ga_thresholds


# ---------------------------------------------------------------------
# Column mappings for HMMER domtblout format
# ---------------------------------------------------------------------
HMMSCAN_COLS = {
    "pfam_full": 1,
    "protein": 3,
    "hmm_len": 2,
    "bitscore": 7,
    "iEvalue": 12,
    "hmm_from": 15,
    "hmm_to": 16,
    "ali_from": 17,
    "ali_to": 18,
}

HMMSEARCH_COLS = {
    "protein": 0,
    "pfam_full": 3,
    "hmm_len": 5,
    "bitscore": 7,
    "iEvalue": 6,
    "hmm_from": 15,
    "hmm_to": 16,
    "ali_from": 17,
    "ali_to": 18,
}


# ---------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------
def parse_pfam_input(pfam_arg: str) -> Set[str]:
    """
    Parse Pfam domain input: a file (one ID per line), a single ID,
    or comma-separated IDs. Returns a set of Pfam IDs (version stripped).
    """
    pfam_ids: Set[str] = set()

    if os.path.isfile(pfam_arg):
        logging.info(f"Reading Pfam IDs from file: {pfam_arg}")
        try:
            with open(pfam_arg, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        pfam_ids.add(line.split(".")[0])
        except IOError as e:
            logging.error(f"Error reading Pfam file {pfam_arg}: {e}")
            sys.exit(1)
    else:
        for pfam_id in pfam_arg.split(","):
            pfam_id = pfam_id.strip()
            if pfam_id:
                pfam_ids.add(pfam_id.split(".")[0])

    if not pfam_ids:
        logging.error("No valid Pfam IDs found in input")
        sys.exit(1)

    return pfam_ids


def extract_genome_name(domtbl_file: str, domtbl_ext: str) -> str:
    """Extract genome name from domtblout filename by removing the extension."""
    basename = os.path.basename(domtbl_file)
    if basename.endswith(domtbl_ext):
        return basename[: -len(domtbl_ext)]
    return basename.rsplit(".", 1)[0]


# ---------------------------------------------------------------------
# File lookup
# ---------------------------------------------------------------------
def create_file_lookup(
    domtbl_dir: str, domtbl_ext: str, genbank_dir: str
) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Build a lookup dict mapping genome basenames to their file paths.
    Keys come from domtblout filenames; GenBank files are matched by prefix.
    """
    file_lookup: Dict[str, Dict[str, Optional[str]]] = {}

    domtbl_files = glob.glob(os.path.join(domtbl_dir, f"*{domtbl_ext}"))

    for domtbl_file in domtbl_files:
        basename = os.path.basename(domtbl_file)
        genome_key = (
            basename[: -len(domtbl_ext)]
            if basename.endswith(domtbl_ext)
            else basename.rsplit(".", 1)[0]
        )
        file_lookup[genome_key] = {"domtbl": domtbl_file, "genbank": None}

    if os.path.exists(genbank_dir):
        gbk_extensions = (".gbk", ".gbff", ".genbank", ".gb")
        for genbank_file in glob.glob(os.path.join(genbank_dir, "*")):
            basename = os.path.basename(genbank_file)
            for genome_key in file_lookup:
                if basename.startswith(genome_key) and any(
                    basename.endswith(ext) for ext in gbk_extensions
                ):
                    file_lookup[genome_key]["genbank"] = genbank_file
                    break

    return file_lookup


# ---------------------------------------------------------------------
# Domtblout parsing & QC
# ---------------------------------------------------------------------
def parse_domtbl(
    domtbl_file: str, num_cols: int, mode: str, domtbl_ext: str
) -> List[Dict]:
    """Parse one domtblout file and return list of hit dictionaries."""
    hits: List[Dict] = []
    genome = extract_genome_name(domtbl_file, domtbl_ext)

    cols_map = HMMSCAN_COLS if mode == "hmmscan" else HMMSEARCH_COLS

    if not os.path.exists(domtbl_file):
        logging.error(f"File not found: {domtbl_file}")
        return hits

    try:
        with open(domtbl_file, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue

                cols = line.rstrip().split()
                if len(cols) < num_cols:
                    continue

                try:
                    pfam_full = cols[cols_map["pfam_full"]]
                    pfam_id = pfam_full.split(".")[0] if "." in pfam_full else pfam_full

                    rec = {
                        "genome": genome,
                        "protein": cols[cols_map["protein"]],
                        "pfam": pfam_id,
                        "bitscore": float(cols[cols_map["bitscore"]]),
                        "iEvalue": float(cols[cols_map["iEvalue"]]),
                        "hmm_from": int(cols[cols_map["hmm_from"]]),
                        "hmm_to": int(cols[cols_map["hmm_to"]]),
                        "ali_from": int(cols[cols_map["ali_from"]]),
                        "ali_to": int(cols[cols_map["ali_to"]]),
                        "hmm_len": int(cols[cols_map["hmm_len"]]),
                    }

                    if (
                        rec["hmm_from"] > rec["hmm_to"]
                        or rec["ali_from"] > rec["ali_to"]
                        or rec["hmm_from"] < 1
                        or rec["ali_from"] < 1
                    ):
                        continue

                    hits.append(rec)
                except (IndexError, ValueError):
                    continue

    except IOError as e:
        logging.error(f"Cannot read file {domtbl_file}: {e}")
        return []

    return hits


def qc_filter(
    hits: List[Dict],
    evalue_thresh: float,
    cov_thresh: float,
    ga_thresholds: Optional[Dict[str, Dict[str, float]]] = None,
    use_ga: bool = False,
) -> List[Dict]:
    """Apply E-value + HMM coverage filter, optionally using GA thresholds."""
    if not hits:
        return []

    passed = []
    ga_used = 0
    fallback_used = 0

    for h in hits:
        try:
            hmm_span = h["hmm_to"] - h["hmm_from"] + 1
            cov = hmm_span / h["hmm_len"] if h["hmm_len"] > 0 else 0

            if cov < cov_thresh:
                continue

            pass_score = False

            if use_ga and ga_thresholds and h["pfam"] in ga_thresholds:
                domain_ga = ga_thresholds[h["pfam"]]["domain_ga"]
                if h["bitscore"] >= domain_ga:
                    pass_score = True
                    ga_used += 1
            else:
                if h["iEvalue"] <= evalue_thresh:
                    pass_score = True
                    fallback_used += 1

            if pass_score:
                h["coverage"] = cov
                passed.append(h)

        except (KeyError, ZeroDivisionError) as e:
            logging.warning(f"QC filter error for hit: {e}")
            continue

    if use_ga:
        logging.info(
            f"QC filter: kept {len(passed):,}/{len(hits):,} hits "
            f"({ga_used:,} via GA, {fallback_used:,} via E-value)"
        )
    else:
        logging.info(f"QC filter: kept {len(passed):,}/{len(hits):,} hits")

    return passed


def resolve_overlaps(
    hits: List[Dict], min_len: int = 15, min_frac: float = 0.25
) -> List[Dict]:
    """Remove substantially overlapping domains on the same protein."""
    if not hits:
        return []

    by_prot: Dict[Tuple[str, str], List[Dict]] = defaultdict(list)
    for h in hits:
        by_prot[(h["genome"], h["protein"])].append(h)

    retained = []
    for group in by_prot.values():
        group.sort(
            key=lambda x: (-x["bitscore"], x["iEvalue"], -(x["ali_to"] - x["ali_from"]))
        )

        accepted = []
        for h in group:
            clash = False
            h_len = h["ali_to"] - h["ali_from"] + 1

            for g in accepted:
                overlap_start = max(h["ali_from"], g["ali_from"])
                overlap_end = min(h["ali_to"], g["ali_to"])
                overlap_len = max(0, overlap_end - overlap_start + 1)

                if overlap_len <= min_len:
                    continue

                g_len = g["ali_to"] - g["ali_from"] + 1
                shorter_len = min(h_len, g_len)

                if shorter_len > 0 and overlap_len / shorter_len >= min_frac:
                    clash = True
                    break

            if not clash:
                accepted.append(h)

        retained.extend(accepted)

    return retained


# ---------------------------------------------------------------------
# Domtblout processing worker (for parallel execution)
# ---------------------------------------------------------------------
def process_domtbl_file(
    domtbl_file: str,
    num_cols: int,
    mode: str,
    evalue: float,
    coverage: float,
    no_overlaps: bool,
    ov_min_len: int,
    ov_min_frac: float,
    domtbl_ext: str,
    ga_thresholds: Optional[Dict[str, Dict[str, float]]] = None,
    use_ga: bool = False,
) -> List[Dict]:
    """Process a single domtblout file through the full pipeline."""
    try:
        hits = parse_domtbl(domtbl_file, num_cols, mode, domtbl_ext)
        if not hits:
            return []

        hits = qc_filter(hits, evalue, coverage, ga_thresholds, use_ga)
        if not hits:
            return []

        if no_overlaps:
            hits = resolve_overlaps(hits, ov_min_len, ov_min_frac)

        return hits
    except Exception as e:
        logging.error(f"Error processing domtblout file {domtbl_file}: {e}")
        return []


# ---------------------------------------------------------------------
# Neighbourhood extraction
# ---------------------------------------------------------------------
def extract_protein_id(feature) -> Optional[str]:
    """Extract protein ID from a CDS feature (locus_tag first for HMMER compat)."""
    for qualifier in ["locus_tag", "protein_id", "gene", "old_locus_tag"]:
        if qualifier in feature.qualifiers:
            return feature.qualifiers[qualifier][0]
    return None


def process_genome_neighbourhoods(
    genome_data: Tuple[str, List[Dict]],
    upstream: int,
    downstream: int,
    file_lookup: Optional[Dict] = None,
) -> Tuple[str, List[Dict]]:
    """Process neighbourhoods for a single genome."""
    genome_name, genome_hits = genome_data

    try:
        gbk_file = None
        if file_lookup and genome_name in file_lookup:
            gbk_file = file_lookup[genome_name].get("genbank")

        if not gbk_file:
            logging.warning(f"GenBank file not found for genome: {genome_name}")
            return genome_name, []

        neighbourhoods = []
        records = list(SeqIO.parse(gbk_file, "genbank"))

        for record in records:
            cds_features = []
            for i, feature in enumerate(record.features):
                if feature.type == "CDS":
                    protein_id = extract_protein_id(feature)
                    if protein_id:
                        product = feature.qualifiers.get("product", [""])[0]
                        notes = feature.qualifiers.get("note", [])
                        cds_features.append(
                            {
                                "index": i,
                                "feature": feature,
                                "protein_id": protein_id,
                                "start": int(feature.location.start),
                                "end": int(feature.location.end),
                                "strand": feature.location.strand,
                                "product": product,
                                "note": "; ".join(notes) if notes else "",
                            }
                        )

            cds_features.sort(key=lambda x: x["start"])
            protein_to_idx = {
                cds["protein_id"]: idx for idx, cds in enumerate(cds_features)
            }

            for hit in genome_hits:
                protein_id = hit["protein"]

                if protein_id not in protein_to_idx:
                    continue

                target_idx = protein_to_idx[protein_id]
                start_idx = max(0, target_idx - upstream)
                end_idx = min(len(cds_features), target_idx + downstream + 1)
                neighbourhood_genes = cds_features[start_idx:end_idx]

                actual_upstream = target_idx - start_idx
                actual_downstream = (end_idx - 1) - target_idx

                neighbourhoods.append(
                    {
                        "genome": genome_name,
                        "record_id": record.id,
                        "record": record,
                        "target_pfam": hit["pfam"],
                        "target_protein": protein_id,
                        "target_gene_idx": target_idx,
                        "neighbourhood_start_idx": start_idx,
                        "neighbourhood_end_idx": end_idx - 1,
                        "genes": neighbourhood_genes,
                        "hit_info": hit,
                        "actual_upstream": actual_upstream,
                        "actual_downstream": actual_downstream,
                        "requested_upstream": upstream,
                        "requested_downstream": downstream,
                    }
                )

        return genome_name, neighbourhoods

    except Exception as e:
        logging.error(f"Error processing genome {genome_name}: {e}")
        return genome_name, []


# ---------------------------------------------------------------------
# GenBank output
# ---------------------------------------------------------------------
def create_neighbourhood_genbank(neighbourhood_info: Dict, output_dir: str) -> str:
    """Create a new GenBank file containing only the neighbourhood genes."""
    genome = neighbourhood_info["genome"]
    pfam = neighbourhood_info["target_pfam"]
    protein = neighbourhood_info["target_protein"]
    record = neighbourhood_info["record"]
    genes = neighbourhood_info["genes"]
    actual_upstream = neighbourhood_info["actual_upstream"]
    actual_downstream = neighbourhood_info["actual_downstream"]

    pfam_output_dir = os.path.join(output_dir, pfam)
    os.makedirs(pfam_output_dir, exist_ok=True)

    safe_protein = re.sub(r"[^\w\-_.]", "_", protein)
    filename = f"{genome}_{pfam}_{safe_protein}_neighbourhood.gbk"
    filepath = os.path.join(pfam_output_dir, filename)

    start_pos = genes[0]["start"]
    end_pos = genes[-1]["end"]
    neighbourhood_seq = record.seq[start_pos:end_pos]

    new_record = SeqRecord(
        neighbourhood_seq,
        id=f"{record.id}_neighbourhood_{start_pos}_{end_pos}",
        description=(
            f"Genomic neighbourhood around {pfam} domain in {protein} from {genome} "
            f"({actual_upstream}+1+{actual_downstream} genes, "
            f"contig: {record.id}, pos: {start_pos}-{end_pos})"
        ),
    )

    new_record.annotations = record.annotations.copy()
    new_record.annotations.setdefault("molecule_type", "DNA")
    new_record.annotations.setdefault("data_file_division", "BCT")
    new_record.annotations.setdefault("topology", "linear")
    if "date" not in new_record.annotations:
        from datetime import datetime

        new_record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()

    for i, gene in enumerate(genes):
        feature = gene["feature"]
        new_start = gene["start"] - start_pos
        new_end = gene["end"] - start_pos
        new_location = FeatureLocation(new_start, new_end, strand=feature.location.strand)

        new_feature = SeqFeature(
            location=new_location,
            type=feature.type,
            qualifiers=feature.qualifiers.copy(),
        )

        notes = new_feature.qualifiers.get("note", []).copy()

        if gene["protein_id"] == protein:
            notes.append(f"TARGET: Contains {pfam} domain")
            notes.append(f"Position in neighbourhood: {actual_upstream + 1} of {len(genes)}")
        else:
            relative_pos = i - actual_upstream
            if relative_pos < 0:
                notes.append(f"Upstream gene: {abs(relative_pos)} genes before target")
            else:
                notes.append(f"Downstream gene: {relative_pos} genes after target")

        new_feature.qualifiers["note"] = notes
        new_record.features.append(new_feature)

    try:
        SeqIO.write(new_record, filepath, "genbank")
    except Exception as e:
        logging.error(f"Failed to write GenBank file {filepath}: {e}")
        raise

    return filepath


# ---------------------------------------------------------------------
# FAA output
# ---------------------------------------------------------------------
def create_pfam_faa_files(
    neighbourhoods: Dict[str, List[Dict]], output_dir: str
) -> None:
    """Create .faa files with protein sequences for each Pfam domain."""
    logging.info("Creating .faa files for each Pfam domain...")

    faa_dir = os.path.join(output_dir, "faa_files")
    os.makedirs(faa_dir, exist_ok=True)

    pfam_genes: Dict[str, Dict[str, Dict]] = defaultdict(dict)

    for genome_neighbourhoods in neighbourhoods.values():
        for neighbourhood in genome_neighbourhoods:
            pfam_domain = neighbourhood["target_pfam"]

            for gene in neighbourhood["genes"]:
                protein_id = gene["protein_id"]
                feature = gene["feature"]

                if "translation" not in feature.qualifiers:
                    continue
                if protein_id in pfam_genes[pfam_domain]:
                    continue

                protein_seq = feature.qualifiers["translation"][0]
                product = gene.get("product", "")

                header = f">{protein_id}"
                if product:
                    header += f" {product}"

                pfam_genes[pfam_domain][protein_id] = {
                    "header": header,
                    "sequence": protein_seq,
                }

    files_created = 0
    total_sequences = 0

    for pfam_domain, genes in pfam_genes.items():
        if not genes:
            continue

        faa_file = os.path.join(faa_dir, f"{pfam_domain}.faa")

        try:
            with open(faa_file, "w") as f:
                for gene in genes.values():
                    f.write(f"{gene['header']}\n")
                    seq = gene["sequence"]
                    for i in range(0, len(seq), 80):
                        f.write(f"{seq[i:i + 80]}\n")

            files_created += 1
            total_sequences += len(genes)
            logging.info(f"Created {faa_file} with {len(genes)} sequences")
        except Exception as e:
            logging.error(f"Error writing .faa file {faa_file}: {e}")

    logging.info(
        f"Created {files_created} .faa files with {total_sequences} total sequences"
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Extract genomic neighbourhoods around Pfam domains from GenBank files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required
    parser.add_argument(
        "--domtbls", required=True, help="Directory containing domtblout files"
    )
    parser.add_argument(
        "--genbanks", required=True, help="Directory containing GenBank files"
    )
    parser.add_argument(
        "--pfam",
        required=True,
        help="Target Pfam domain ID(s): single ID (PF00001), "
        "comma-separated (PF00001,PF00002), or path to file with one ID per line",
    )

    # Output
    parser.add_argument(
        "--output", default="neighbourhoods", help="Output directory"
    )

    # File extensions
    parser.add_argument(
        "--domtbl-ext", default=".domtblout", help="File extension for domtblout files"
    )

    # Neighbourhood size
    parser.add_argument(
        "--upstream", type=int, default=5, help="Number of genes upstream to include"
    )
    parser.add_argument(
        "--downstream", type=int, default=5, help="Number of genes downstream to include"
    )

    # HMMER mode
    parser.add_argument(
        "--mode",
        choices=["hmmscan", "hmmsearch"],
        default="hmmscan",
        help="HMMER mode: hmmscan (seqs vs HMMs) or hmmsearch (HMMs vs seqs)",
    )
    parser.add_argument(
        "--num-cols", type=int, default=23, help="Expected columns in domtblout files"
    )

    # QC thresholds
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-3,
        help="Maximum independent E-value threshold (fallback when GA not available)",
    )
    parser.add_argument(
        "--coverage", type=float, default=0.5, help="Minimum HMM coverage fraction"
    )

    # Overlap resolution
    parser.add_argument(
        "--no-overlaps",
        action="store_true",
        help="Apply overlap resolution to remove conflicting domains",
    )
    parser.add_argument(
        "--ov-min-len", type=int, default=15, help="Minimum overlap length"
    )
    parser.add_argument(
        "--ov-min-frac", type=float, default=0.25, help="Minimum overlap fraction"
    )

    # Neighbourhood filter
    parser.add_argument(
        "--min-genes",
        type=int,
        default=0,
        help="Minimum total genes required in neighbourhood (0 = no minimum)",
    )

    # Parallel processing
    parser.add_argument(
        "--cpu", type=int, default=1, help="Number of CPU cores for parallel processing"
    )

    # GA thresholds
    parser.add_argument(
        "--ga-file", type=str, default=None, help="Path to Pfam-A.hmm file for GA thresholds"
    )
    parser.add_argument(
        "--use-ga",
        action="store_true",
        help="Use GA thresholds when available (requires --ga-file)",
    )

    # Optional outputs
    parser.add_argument(
        "--create-faa",
        action="store_true",
        help="Create .faa files with protein sequences for each Pfam domain",
    )

    # Verbosity
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity: -v = INFO, -vv = DEBUG",
    )

    args = parser.parse_args()

    # Logging
    if args.verbose >= 2:
        level = logging.DEBUG
    elif args.verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Validate
    if args.cpu < 1:
        logging.error("CPU count must be at least 1")
        sys.exit(1)

    if args.use_ga and not args.ga_file:
        logging.error("--use-ga requires --ga-file")
        sys.exit(1)

    system_cpus = multiprocessing.cpu_count()
    if args.cpu > system_cpus:
        logging.warning(
            f"Requested {args.cpu} CPUs but system has {system_cpus}. "
            f"Using {system_cpus}."
        )
        args.cpu = system_cpus

    if not os.path.isdir(args.domtbls):
        logging.error(f"Domtblout directory not found: {args.domtbls}")
        sys.exit(1)

    if not os.path.isdir(args.genbanks):
        logging.error(f"GenBank directory not found: {args.genbanks}")
        sys.exit(1)

    # GA thresholds
    ga_thresholds = None
    if args.use_ga:
        ga_thresholds = parse_ga_thresholds(args.ga_file)
        if not ga_thresholds:
            logging.error("No GA thresholds found - check your --ga-file path")
            sys.exit(1)

    # Parse Pfam input
    target_pfams = parse_pfam_input(args.pfam)
    logging.info(f"Target Pfam domains: {sorted(target_pfams)}")

    os.makedirs(args.output, exist_ok=True)

    # Build file lookup
    logging.info("Building file lookup...")
    file_lookup = create_file_lookup(args.domtbls, args.domtbl_ext, args.genbanks)
    genomes_with_gbk = sum(1 for v in file_lookup.values() if v.get("genbank"))
    logging.info(
        f"File lookup: {len(file_lookup)} genomes, "
        f"{genomes_with_gbk} with GenBank files"
    )

    # Find domtblout files
    pattern = os.path.join(args.domtbls, f"*{args.domtbl_ext}")
    domtbl_files = sorted(glob.glob(pattern))

    if not domtbl_files:
        logging.error(f"No domtblout files found matching {pattern}")
        sys.exit(1)

    logging.info(f"Found {len(domtbl_files)} domtblout files")

    # ---- Parse domtblout files ----
    domtbl_worker = partial(
        process_domtbl_file,
        num_cols=args.num_cols,
        mode=args.mode,
        evalue=args.evalue,
        coverage=args.coverage,
        no_overlaps=args.no_overlaps,
        ov_min_len=args.ov_min_len,
        ov_min_frac=args.ov_min_frac,
        domtbl_ext=args.domtbl_ext,
        ga_thresholds=ga_thresholds,
        use_ga=args.use_ga,
    )

    all_hits: List[Dict] = []

    if args.cpu > 1:
        with ProcessPoolExecutor(max_workers=args.cpu) as executor:
            futures = {executor.submit(domtbl_worker, f): f for f in domtbl_files}
            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Processing domtblout files",
                unit="file",
            ):
                all_hits.extend(future.result())
    else:
        for domtbl_file in tqdm(
            domtbl_files, desc="Processing domtblout files", unit="file"
        ):
            all_hits.extend(domtbl_worker(domtbl_file))

    if not all_hits:
        logging.error("No hits found after filtering")
        sys.exit(1)

    logging.info(f"Total hits after QC: {len(all_hits)}")

    # ---- Filter for target Pfams ----
    target_hits = [h for h in all_hits if h["pfam"] in target_pfams]

    if not target_hits:
        logging.error(f"No hits found for target Pfam domains: {sorted(target_pfams)}")
        sys.exit(1)

    hits_per_pfam: Dict[str, int] = defaultdict(int)
    for h in target_hits:
        hits_per_pfam[h["pfam"]] += 1

    logging.info(f"Found {len(target_hits)} hits for target Pfam domains:")
    for pfam in sorted(target_pfams):
        logging.info(f"  {pfam}: {hits_per_pfam.get(pfam, 0)} hits")

    # ---- Extract neighbourhoods ----
    hits_by_genome: Dict[str, List[Dict]] = defaultdict(list)
    for h in target_hits:
        hits_by_genome[h["genome"]].append(h)

    genome_data_list = list(hits_by_genome.items())

    worker_func = partial(
        process_genome_neighbourhoods,
        upstream=args.upstream,
        downstream=args.downstream,
        file_lookup=file_lookup,
    )

    neighbourhoods: Dict[str, List[Dict]] = {}

    if args.cpu > 1:
        with ProcessPoolExecutor(max_workers=args.cpu) as executor:
            futures = {
                executor.submit(worker_func, gd): gd[0] for gd in genome_data_list
            }
            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Extracting neighbourhoods",
                unit="genome",
            ):
                genome_name, genome_nbhs = future.result()
                if genome_nbhs:
                    neighbourhoods[genome_name] = genome_nbhs
    else:
        for genome_data in tqdm(
            genome_data_list, desc="Extracting neighbourhoods", unit="genome"
        ):
            genome_name, genome_nbhs = worker_func(genome_data)
            if genome_nbhs:
                neighbourhoods[genome_name] = genome_nbhs

    if not neighbourhoods:
        logging.error("No genomic neighbourhoods found")
        sys.exit(1)

    total_found = sum(len(nbhs) for nbhs in neighbourhoods.values())
    logging.info(f"Found {total_found} neighbourhoods across {len(neighbourhoods)} genomes")

    # ---- Apply min-genes filter ----
    if args.min_genes > 0:
        total_before = total_found
        filtered = {}
        for gn, gnbhs in neighbourhoods.items():
            kept = [n for n in gnbhs if len(n["genes"]) >= args.min_genes]
            if kept:
                filtered[gn] = kept
        neighbourhoods = filtered
        total_after = sum(len(nbhs) for nbhs in neighbourhoods.values())

        logging.info(
            f"Filtered neighbourhoods: {total_after}/{total_before} "
            f"have >= {args.min_genes} genes"
        )
        if not neighbourhoods:
            logging.error(f"No neighbourhoods have >= {args.min_genes} genes")
            sys.exit(1)

    # ---- Write GenBank output ----
    total_neighbourhoods = sum(len(nbhs) for nbhs in neighbourhoods.values())
    logging.info(f"Creating {total_neighbourhoods} neighbourhood GenBank files...")

    all_nbhs = [n for nbhs in neighbourhoods.values() for n in nbhs]

    created_files: List[str] = []
    failed_files = 0
    files_by_pfam: Dict[str, List[str]] = defaultdict(list)

    for neighbourhood in tqdm(all_nbhs, desc="Writing GenBank files", unit="file"):
        try:
            filepath = create_neighbourhood_genbank(neighbourhood, args.output)
            created_files.append(filepath)
            files_by_pfam[neighbourhood["target_pfam"]].append(filepath)
        except Exception as e:
            failed_files += 1
            logging.error(
                f"Error creating file for "
                f"{neighbourhood.get('genome')}:{neighbourhood.get('target_pfam')}:"
                f"{neighbourhood.get('target_protein')} - {e}"
            )

    # ---- FAA output ----
    if args.create_faa:
        try:
            create_pfam_faa_files(neighbourhoods, args.output)
        except Exception as e:
            logging.error(f"Error creating .faa files: {e}")

    # ---- Summary ----
    genes_per_nbh = [len(n["genes"]) for n in all_nbhs]
    avg_genes = sum(genes_per_nbh) / len(genes_per_nbh) if genes_per_nbh else 0

    print(f"\n{'=' * 60}")
    print("PFAM NEIGHBOURHOOD EXTRACTION SUMMARY")
    print(f"{'=' * 60}")
    print(f"Target Pfam domain(s): {', '.join(sorted(target_pfams))}")
    print(f"Processed domtblout files: {len(domtbl_files)}")
    print(f"Total hits after filtering: {len(all_hits)}")
    print(f"Hits for target Pfam(s): {len(target_hits)}")
    print(f"Genomes with target hits: {len(hits_by_genome)}")
    print(f"Neighbourhoods found: {total_found}")
    if args.min_genes > 0:
        print(f"Neighbourhoods after min-genes filter: {total_neighbourhoods}")
    print(f"GenBank files created: {len(created_files)}")
    if failed_files > 0:
        print(f"Failed file creations: {failed_files}")
    print(f"Output directory: {os.path.abspath(args.output)}")
    print(f"Neighbourhood size: {args.upstream}+1+{args.downstream} genes")
    if genes_per_nbh:
        print(f"Average genes per neighbourhood: {avg_genes:.1f}")

    print(f"\nFiles per Pfam domain:")
    for pfam in sorted(files_by_pfam):
        print(f"  {pfam}: {len(files_by_pfam[pfam])} files")

    if args.create_faa:
        print(f"\n.faa files: {os.path.join(args.output, 'faa_files')}")

    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
