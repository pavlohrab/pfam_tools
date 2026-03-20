#!/usr/bin/env python3
"""
pfam_mapping.py – Pfam → description + GO + GO-Slim annotation.

Modes
-----
**Single list:**    annotate one set of Pfam IDs.
**Feature sets:**   annotate several sets, then run universality analysis.

CLI
---
    # single list
    python pfam_mapping.py -i pfams.txt -o results/ -d db/

    # multiple feature sets
    python pfam_mapping.py --feature-sets sets.tsv -o results/ -d db/
"""

from __future__ import annotations
import argparse
import gzip
import pathlib
import re
import time
import urllib.request
from typing import Dict, List
from collections import defaultdict
import numpy as np

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import requests
from goatools.mapslim import mapslim
from goatools.obo_parser import GODag
from pfam2go import pfam2go
import pfam2go as _pfam2go_mod
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Monkey-patch pfam2go: original URL now returns 403.
# ---------------------------------------------------------------------------
_orig_raw_data_to_frame = _pfam2go_mod._raw_data_to_frame

def _patched_data_init():
    import urllib.request as _ur
    urls = [
        "https://raw.githubusercontent.com/geneontology/go-site/master/metadata/external2go/pfam2go",
        "http://current.geneontology.org/ontology/external2go/pfam2go",
    ]
    for url in urls:
        try:
            req = _ur.Request(url, headers={"User-Agent": "pfam-tools/1.0"})
            raw = _ur.urlopen(req, timeout=30).read().decode("utf-8")
            return _orig_raw_data_to_frame(raw)
        except Exception:
            continue
    raise RuntimeError("Could not download pfam2go mapping from any known URL")

_pfam2go_mod._data_init = _patched_data_init

# ---------------------------------------------------------------------------
# Constants & session
# ---------------------------------------------------------------------------
PFAM_TSV_FTP = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
GO_BASIC_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GOSLIM_URL   = "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"

session = requests.Session()
session.headers.update({"User-Agent": "pfam-tools/1.0"})


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Annotate Pfam IDs with GO / GO-Slim.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-i", "--input", default="pfams.txt",
                    help="Pfam IDs: plain-text list OR pfam_presence matrix TSV "
                         "from pfam_scan.py (auto-detected; adds count column)")
    ap.add_argument("-o", "--outdir", default=".",
                    help="Output directory")
    ap.add_argument("-d", "--database", default=".",
                    help="Cache directory for OBO / TSV database files")
    ap.add_argument("--pfam-tsv", dest="pfam_tsv", default=None,
                    help="Path to Pfam-A.clans.tsv[.gz] (auto-downloaded if missing)")
    ap.add_argument("--feature-sets", dest="feature_sets", default=None,
                    help="TSV file with columns: set_name, features "
                         "(comma-separated Pfam IDs per row)")
    return ap.parse_args()


def download_if_missing(url: str, dest: pathlib.Path) -> pathlib.Path:
    if dest.exists():
        return dest
    print(f"Downloading {dest.name} …")
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "pfam-tools/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as fh:
        fh.write(resp.read())
    return dest


def load_pfam_descriptions(tsv_path: pathlib.Path) -> Dict[str, str]:
    opener = gzip.open if tsv_path.suffix == ".gz" else open
    mapping: Dict[str, str] = {}
    with opener(tsv_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                acc, *_, descr = parts[:5]
                mapping[acc] = descr
    return mapping


def pfam_desc_api(pfam_id: str, tries: int = 4) -> str:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}/?page_size=1"
    for n in range(tries):
        try:
            r = session.get(url, timeout=(3, 30))
            if r.status_code == 408:
                time.sleep(int(r.headers.get("Retry-After", 10)))
                continue
            if r.ok:
                js = r.json()
                if js.get("results"):
                    md = js["results"][0].get("metadata", {})
                    return md.get("name") or (md.get("description", ["NA"])[0][:120])
            if r.status_code >= 500:
                time.sleep(2 ** n)
        except (requests.Timeout, requests.ConnectionError):
            time.sleep(2 ** n)
    return "NA"


def flatten_any(x) -> List[str]:
    out: List[str] = []
    stack = [x]
    while stack:
        item = stack.pop()
        if isinstance(item, (set, list, tuple, frozenset)):
            stack.extend(item)
        elif pd.isna(item) or item in ("", "set()"):
            continue
        else:
            out.append(str(item))
    return out


def slim_id_list(go_id: str, full: GODag, slim: GODag) -> List[str]:
    raw = mapslim(go_id, full, slim)
    return sorted(set(flatten_any(raw)))


def _load_input(path: str):
    """
    Auto-detect input format: pfam_presence matrix (TSV) or plain Pfam list.

    Returns:
        (pfams, pfam_counts) where pfam_counts is a dict {PFxxxxx: int} or None.
    """
    inp = pathlib.Path(path)

    # Peek at the first few lines to detect format (avoid reading entire file)
    with open(inp) as fh:
        head = [fh.readline() for _ in range(12)]
    head = [l for l in head if l.strip()]

    if len(head) >= 2 and "\t" in head[0]:
        first_cols = [l.split("\t", 1)[0] for l in head[1:]]
        pfam_like = sum(1 for c in first_cols if re.match(r"PF\d{5}", c, re.I))
        if pfam_like >= len(first_cols) * 0.8:
            print(f"Detected pfam_presence matrix in {inp.name}")
            # Read index column, then load numeric data with numpy
            pfam_ids = []
            with open(inp) as fh:
                header = fh.readline()
                n_samples = header.count("\t")
                for line in fh:
                    pfam_ids.append(line.split("\t", 1)[0])
            ncols = list(range(1, n_samples + 1))
            mat = np.loadtxt(inp, delimiter="\t", skiprows=1,
                             usecols=ncols, dtype=np.int8)
            row_sums = mat.sum(axis=1)
            pfam_counts = dict(zip(pfam_ids, row_sums.tolist()))
            pfams = sorted(pfam_counts.keys())
            print(f"  {len(pfams)} Pfam families, "
                  f"{n_samples} samples, "
                  f"total count {sum(pfam_counts.values()):,}")
            return pfams, pfam_counts

    # Fallback: plain text with Pfam IDs
    raw = inp.read_text()
    pfams = sorted(set(re.findall(r"(?:PFAM_)?(PF\d{5})", raw, flags=re.I)))
    return pfams, None


def load_feature_sets(path: str) -> Dict[str, List[str]]:
    df = pd.read_csv(path, sep="\t")
    result: Dict[str, List[str]] = {}
    for _, row in df.iterrows():
        name = row["set_name"]
        pfams: List[str] = []
        for tok in row["features"].split(","):
            pfams.extend(re.findall(r"(?:PFAM_)?(PF\d{5})", tok, flags=re.I))
        result[name] = sorted(set(pfams))
    return result


def generate_universality_analysis(base_dir: pathlib.Path) -> None:
    feature_sets = [
        d.name for d in base_dir.iterdir()
        if d.is_dir() and (d / "pfam_go_annotation.tsv").exists()
    ]
    if len(feature_sets) < 2:
        print("Need at least 2 feature sets for universality analysis")
        return

    print(f"Generating universality analysis for {len(feature_sets)} feature sets")

    go_terms: Dict[str, set] = defaultdict(set)
    goslim_terms: Dict[str, set] = defaultdict(set)
    go_names: Dict[str, str] = {}

    for fs in feature_sets:
        try:
            df = pd.read_csv(base_dir / fs / "pfam_go_annotation.tsv", sep="\t")
        except Exception as e:
            print(f"Error processing {fs}: {e}")
            continue
        for _, row in df.iterrows():
            ga, gn, gs = row.get("go_accession"), row.get("go_name"), row.get("go_slim_name")
            if pd.notna(ga) and ga:
                go_terms[ga].add(fs)
                go_names[ga] = gn
            if pd.notna(gs) and gs:
                goslim_terms[gs].add(fs)

    n_sets = len(feature_sets)
    if go_terms:
        rows = sorted(
            [{"go_accession": ga, "go_name": go_names.get(ga, ""),
              "universality_count": len(fss), "total_feature_sets": n_sets,
              "universality_pct": round(len(fss) / n_sets * 100, 1),
              "feature_sets": ";".join(sorted(fss))}
             for ga, fss in go_terms.items()],
            key=lambda r: -r["universality_count"],
        )
        pd.DataFrame(rows).to_csv(base_dir / "go_terms_universality.tsv", sep="\t", index=False)
        print(f"GO terms universality → {base_dir / 'go_terms_universality.tsv'}")

    if goslim_terms:
        rows = sorted(
            [{"goslim_term": t, "universality_count": len(fss),
              "total_feature_sets": n_sets,
              "universality_pct": round(len(fss) / n_sets * 100, 1),
              "feature_sets": ";".join(sorted(fss))}
             for t, fss in goslim_terms.items()],
            key=lambda r: -r["universality_count"],
        )
        pd.DataFrame(rows).to_csv(base_dir / "goslim_terms_universality.tsv", sep="\t", index=False)
        print(f"GOslim terms universality → {base_dir / 'goslim_terms_universality.tsv'}")

    print(f"Summary: {len(go_terms)} GO terms, {len(goslim_terms)} GOslim terms analyzed")


def process_feature_set(pfams: List[str], outdir: pathlib.Path, desc_map: Dict[str, str],
                        godag_full: GODag, godag_slim: GODag,
                        pfam_counts: Dict[str, int] | None = None) -> None:
    if not pfams:
        print("No Pfam IDs found in feature set.")
        return

    print(f"Processing {len(pfams)} Pfams")

    pf_desc = {pf: desc_map.get(pf, "NA") for pf in pfams}
    missing = [pf for pf, d in pf_desc.items() if d == "NA"]
    if missing:
        print(f"Fetching {len(missing)} missing descriptions via API")
        for pf in tqdm(missing):
            pf_desc[pf] = pfam_desc_api(pf)

    rows = []
    for pf in pfams:
        row = {"pfam_id": pf, "pfam_desc": pf_desc[pf]}
        if pfam_counts is not None:
            row["count"] = pfam_counts.get(pf, 0)
        rows.append(row)
    pd.DataFrame(rows).to_csv(outdir / "pfam_descriptions_only.tsv", sep="\t", index=False)
    print(f"Pfam descriptions (all) → {outdir / 'pfam_descriptions_only.tsv'}")

    df = pfam2go(pfams)
    if df.empty:
        print("No GO annotations found for any Pfams in this set.")
        return

    df.columns = df.columns.str.lower().str.replace(r"[^\w]+", "_", regex=True).str.strip("_")
    pfcol = next(c for c in df.columns if c.startswith("pfam"))
    gocol = next(c for c in df.columns if "go_accession" in c)

    df["go_slim_ids"] = df[gocol].apply(lambda g: slim_id_list(g, godag_full, godag_slim))
    df["go_slim_name"] = df["go_slim_ids"].apply(
        lambda ids: "; ".join(
            godag_slim[i].name if i in godag_slim else godag_full.get(i, ("",)).name
            for i in ids
        )
    )
    df = df.assign(
        pfam_desc=lambda d: d[pfcol].map(pf_desc),
        go_name=lambda d: d[gocol].map(lambda g: godag_full[g].name if g in godag_full else "NA"),
    )

    out1 = outdir / "pfam_go_annotation.tsv"
    out_cols = [pfcol, "pfam_desc", gocol, "go_name", "go_slim_name"]
    if pfam_counts is not None:
        df["count"] = df[pfcol].map(pfam_counts).fillna(0).astype(int)
        out_cols.insert(2, "count")
    df[out_cols].to_csv(out1, sep="\t", index=False)
    print(f"Annotation table (with GO) → {out1}")

    df_long = (
        df[[pfcol, "go_slim_ids"]]
        .explode("go_slim_ids")
        .rename(columns={"go_slim_ids": "goslim"})
        .drop_duplicates([pfcol, "goslim"])
    )
    summary = df_long.groupby("goslim")[pfcol].nunique().sort_values(ascending=False)
    summary_df = (
        summary.reset_index(name="pfam_count")
        .assign(term_name=lambda d: d["goslim"].map(
            lambda g: godag_slim[g].name if g in godag_slim else godag_full.get(g, ("NA",)).name))
    )[["goslim", "term_name", "pfam_count"]]
    summary_df.to_csv(outdir / "pfam_goslim_summary.tsv", sep="\t", index=False)
    print(f"GO-Slim summary → {outdir / 'pfam_goslim_summary.tsv'}")

    out_text = outdir / "pfam_goslim_explanation.txt"
    with open(out_text, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("GO-SLIM FUNCTIONAL CATEGORIES - EXPLANATION\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total Pfam domains analyzed: {len(pfams)}\n")
        f.write(f"Pfam domains with GO annotations: {df[pfcol].nunique()}\n")
        f.write(f"Unique GO-Slim categories found: {len(summary_df)}\n\n")
        f.write("-" * 80 + "\n\n")
        for _, row in summary_df.iterrows():
            gid = row["goslim"]
            defn = "No definition available"
            for dag in (godag_slim, godag_full):
                if gid in dag and hasattr(dag[gid], "defn"):
                    defn = dag[gid].defn
                    break
            f.write(f"[{row['term_name']}] ({gid})\n")
            f.write(f"  Pfam domains: {row['pfam_count']}\n")
            f.write(f"  Description: {defn}\n\n")
    print(f"GO-Slim explanations → {out_text}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [3, 1]})
    summary_df.set_index("term_name")["pfam_count"].sort_values(ascending=False).plot(
        kind="bar", ax=ax1, color="steelblue")
    ax1.set_ylabel("Pfam domains count", fontsize=12)
    ax1.set_xlabel("GO-Slim term", fontsize=12)
    ax1.set_title("GO-Slim functional categories", fontsize=14, fontweight="bold")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right", fontsize=9)
    ax1.grid(axis="y", alpha=0.3)

    annotated = df[pfcol].nunique()
    unannotated = len(pfams) - annotated
    ax2.pie([annotated, unannotated],
            labels=[f"Annotated\n({annotated})", f"Not annotated\n({unannotated})"],
            autopct="%1.1f%%", startangle=90, colors=["#66c2a5", "#fc8d62"],
            textprops={"fontsize": 10})
    ax2.set_title("GO annotation coverage", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(outdir / "pfam_goslim_summary.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot → {outdir / 'pfam_goslim_summary.png'}")


def main():
    args = parse_args()
    outdir = pathlib.Path(args.outdir).expanduser().resolve()
    dbdir = pathlib.Path(args.database).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    dbdir.mkdir(parents=True, exist_ok=True)

    pfam_tsv = pathlib.Path(args.pfam_tsv) if args.pfam_tsv else dbdir / "Pfam-A.clans.tsv.gz"
    pfam_tsv = download_if_missing(PFAM_TSV_FTP, pfam_tsv)
    desc_map = load_pfam_descriptions(pfam_tsv)
    print(f"Loaded {len(desc_map)} Pfam descriptions")

    full_obo = download_if_missing(GO_BASIC_URL, dbdir / "go-basic.obo")
    slim_obo = download_if_missing(GOSLIM_URL, dbdir / "goslim_generic.obo")
    godag_full = GODag(str(full_obo))
    godag_slim = GODag(str(slim_obo))

    if args.feature_sets:
        feature_sets = load_feature_sets(args.feature_sets)
        print(f"Loaded {len(feature_sets)} feature sets")
        for set_name, pfams in feature_sets.items():
            print(f"\n=== Processing feature set: {set_name} ===")
            set_outdir = outdir / set_name
            set_outdir.mkdir(parents=True, exist_ok=True)
            process_feature_set(pfams, set_outdir, desc_map, godag_full, godag_slim)
        print(f"\n=== Generating universality analysis ===")
        generate_universality_analysis(outdir)
    else:
        pfams, pfam_counts = _load_input(args.input)
        if not pfams:
            raise SystemExit("No Pfam IDs found.")
        print(f"{len(pfams)} Pfams to process")
        process_feature_set(pfams, outdir, desc_map, godag_full, godag_slim,
                            pfam_counts=pfam_counts)


if __name__ == "__main__":
    main()
