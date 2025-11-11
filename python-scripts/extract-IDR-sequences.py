#!/usr/bin/env python3

import argparse
import ast
import sys
import pandas as pd


def get_args():
    p = argparse.ArgumentParser(
        description=(
            "Extract IDR sequences from IDR table and write them as FASTA "
            "records keyed by IDR_ID."
        )
    )
    p.add_argument(
        "--input",
        required=True,
        help="CSV file with IDR annotations (one row per IDR_ID).",
    )
    p.add_argument(
        "--output-fasta",
        required=True,
        help="Output FASTA file with one record per IDR_ID.",
    )
    return p.parse_args()


def parse_list_field(value, colname, row_idx):
    """
    Safely parse a stringified Python list from the CSV.
    e.g. "['375-451']" -> ['375-451']
    """
    if pd.isna(value):
        return []
    s = str(value).strip()
    if s == "":
        return []
    try:
        parsed = ast.literal_eval(s)
    except (SyntaxError, ValueError) as e:
        print(
            f"[WARN] Could not parse {colname} in row {row_idx}: {value!r} ({e})",
            file=sys.stderr,
        )
        return []
    if not isinstance(parsed, (list, tuple)):
        print(
            f"[WARN] Parsed {colname} in row {row_idx} is not a list/tuple: {parsed!r}",
            file=sys.stderr,
        )
        return []
    return list(parsed)


def main():
    args = get_args()

    df = pd.read_csv(args.input)

    required_cols = [
        "IDR_ID",
        "IDR_START",
        "IDR_END",
        "idr_ranges",
        "idr_sequences",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: Missing required columns: {missing}")

    n_written = 0

    with open(args.output_fasta, "w") as fasta:
        for idx, row in df.iterrows():
            idr_id = str(row["IDR_ID"])
            try:
                start = int(row["IDR_START"])
                end = int(row["IDR_END"])
            except ValueError:
                print(
                    f"[WARN] Non-integer IDR_START/IDR_END in row {idx} for IDR_ID={idr_id}",
                    file=sys.stderr,
                )
                continue

            target_range = f"{start}-{end}"

            # Parse the list-like columns
            ranges = parse_list_field(row["idr_ranges"], "idr_ranges", idx)
            seqs = parse_list_field(row["idr_sequences"], "idr_sequences", idx)

            if not ranges or not seqs:
                print(
                    f"[WARN] Empty idr_ranges/idr_sequences for IDR_ID={idr_id} (row {idx})",
                    file=sys.stderr,
                )
                continue

            if len(ranges) != len(seqs):
                print(
                    f"[WARN] Length mismatch in row {idx} for IDR_ID={idr_id}: "
                    f"{len(ranges)} ranges vs {len(seqs)} sequences",
                    file=sys.stderr,
                )

            # Find the index where the range matches IDR_START-IDR_END
            try:
                i = ranges.index(target_range)
            except ValueError:
                # No exact match; report and skip
                print(
                    f"[WARN] No matching range {target_range} in idr_ranges for "
                    f"IDR_ID={idr_id} (row {idx}); ranges={ranges}",
                    file=sys.stderr,
                )
                continue

            if i >= len(seqs):
                print(
                    f"[WARN] Index {i} for IDR_ID={idr_id} out of bounds for idr_sequences "
                    f"(len={len(seqs)}); ranges={ranges}",
                    file=sys.stderr,
                )
                continue

            seq = str(seqs[i]).strip()
            if not seq:
                print(
                    f"[WARN] Empty sequence for IDR_ID={idr_id} at range {target_range}",
                    file=sys.stderr,
                )
                continue

            # Write FASTA record
            # You can customize the header to include dataset/Gene if desired.
            fasta.write(f">{idr_id}\n")
            # Wrap sequence to 60 or 80 chars/line if you like; here we write as a single line
            fasta.write(seq + "\n")

            n_written += 1

    print(f"[INFO] Wrote {n_written} IDR sequences to {args.output_fasta}")


if __name__ == "__main__":
    main()
