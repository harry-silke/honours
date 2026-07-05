"""
build_metadata_from_master.py
-----------------------------
Reads the Heat_stress master workbook and writes metadata/crosses.csv,
keyed by Cross ID, for the aggregator to join. The numeric part of each
photo id (e.g. 103 in "103B") is the Cross ID.

Output columns
--------------
    cross_id, parent1_culture, parent2_culture,
    parent1_pop, parent2_pop, cross_type   (within / between)

Population is derived from the parent culture-ID prefix via SITE_TO_POP.
Note: most crosses are BETWEEN populations, so there is no single
population per cross; both parent populations are kept.

Usage
-----
    python build_metadata_from_master.py <workbook.xlsx> --out <root>/metadata
"""

import argparse
import os
import re

import pandas as pd

SITE_TO_POP = {
    "LSI": "Little Spectacle Island",
    "IOC": "Isle of Caves",
    "SHL": "Shelly Point",
    "BBB": "Beerbarrel Beach",
    "ACT": "Actaeon Island",
    "FSH": "Fishers Point",
    "ULV": "Black Jack Rock",
    "DEV": "Devonport breakwater",
    "DIP": "Petrel Point",
    "SDP": "Sandpatch Point reef",
    "PHK": "Stable Bay",
}


def site_of(culture_id):
    m = re.match(r"([A-Za-z]+)", str(culture_id))
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workbook")
    ap.add_argument("--out", default="metadata")
    ap.add_argument("--sheet", default="Crosses")
    args = ap.parse_args()

    cr = pd.ExcelFile(args.workbook).parse(args.sheet)
    cr = cr[cr["Cross ID"].notna()].copy()
    cr["cross_id"] = cr["Cross ID"].astype(int)

    p1 = "Parent 1 - culture ID"
    p2 = "Parent 2 - culture ID"
    out = pd.DataFrame({
        "cross_id": cr["cross_id"],
        "parent1_culture": cr[p1].astype(str),
        "parent2_culture": cr[p2].astype(str),
    })
    out["parent1_pop"] = out["parent1_culture"].map(
        lambda x: SITE_TO_POP.get(site_of(x)))
    out["parent2_pop"] = out["parent2_culture"].map(
        lambda x: SITE_TO_POP.get(site_of(x)))
    out["cross_type"] = [
        "within" if a == b else "between"
        for a, b in zip(out["parent1_pop"], out["parent2_pop"])
    ]

    unmapped = sorted(set(
        site_of(c) for c in pd.concat([out["parent1_culture"],
                                       out["parent2_culture"]])
        if SITE_TO_POP.get(site_of(c)) is None
    ))
    if unmapped:
        print(f"WARNING: unmapped site codes (add to SITE_TO_POP): {unmapped}")

    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, "crosses.csv")
    out.sort_values("cross_id").to_csv(path, index=False)
    print(f"Wrote {len(out)} crosses to {path}")
    print(out["cross_type"].value_counts().to_string())


if __name__ == "__main__":
    main()