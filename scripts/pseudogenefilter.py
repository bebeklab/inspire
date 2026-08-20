import csv
from pathlib import Path

input_path = Path(
    "networks/Blossom_7_31_rulecount_cutoff_v2.csv"
)

output_path = Path(
    "networks/Blossom_7_31_rulecount_cutoff_v2_protein_nodes_only.csv"
)

excluded_path = Path(
    "reports/excluded_paths_nonprotein_nodes.csv"
)

excluded_symbols = {
    "PMS2CL",
}

kept_rows = 0
excluded_rows = 0
excluded_paths = set()

output_path.parent.mkdir(parents=True, exist_ok=True)
excluded_path.parent.mkdir(parents=True, exist_ok=True)

with input_path.open(
    newline="",
    encoding="utf-8-sig",
) as input_handle:
    reader = csv.DictReader(input_handle)

    if reader.fieldnames is None:
        raise RuntimeError("Input CSV has no header")

    required = {
        "gene1",
        "gene2",
        "path",
        "score",
        "petal",
    }

    missing = required - set(reader.fieldnames)

    if missing:
        raise RuntimeError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as output_handle, excluded_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as excluded_handle:

        writer = csv.DictWriter(
            output_handle,
            fieldnames=reader.fieldnames,
        )

        excluded_writer = csv.DictWriter(
            excluded_handle,
            fieldnames=[
                "source_row",
                "reason",
                *reader.fieldnames,
            ],
        )

        writer.writeheader()
        excluded_writer.writeheader()

        for source_row, row in enumerate(
            reader,
            start=2,
        ):
            path_text = row.get("path", "")

            matched_symbols = sorted(
                symbol
                for symbol in excluded_symbols
                if symbol in path_text.split()
            )

            if matched_symbols:
                excluded_rows += 1
                excluded_paths.add(
                    (
                        row.get("petal", ""),
                        path_text,
                    )
                )

                excluded_writer.writerow(
                    {
                        "source_row": source_row,
                        "reason": (
                            "nonprotein_node:"
                            + ",".join(matched_symbols)
                        ),
                        **row,
                    }
                )
                continue

            writer.writerow(row)
            kept_rows += 1

print(f"Input:          {input_path}")
print(f"Filtered:       {output_path}")
print(f"Audit report:   {excluded_path}")
print(f"Kept rows:      {kept_rows}")
print(f"Excluded rows:  {excluded_rows}")
print(f"Excluded paths: {len(excluded_paths)}")

for petal, path_text in sorted(excluded_paths):
    print(f"  {petal}: {path_text}")
