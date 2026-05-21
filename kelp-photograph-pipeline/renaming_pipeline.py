"""
Jar-photo pipeline: ARW -> JPEG + TIFF -> Claude reads jar IDs -> rename TIFFs.

Run as three separate stages. The CSV review between stage 2 and stage 3 is a
HARD STOP: stage 3 refuses to run unless you have looked at the CSV and marked
it approved (see --i-have-reviewed below).

Pass the raws folder as an argument. jpegs/, tiffs/ and rename_map.csv are
created alongside it (in the same parent), so each shoot stays self-contained:

    python renaming_pipeline.py convert /path/to/shoot1/raws
    python renaming_pipeline.py extract /path/to/shoot1/raws
    # ---- now open shoot1/rename_map.csv and check it by hand ----
    python renaming_pipeline.py apply /path/to/shoot1/raws --i-have-reviewed

If you omit the path it falls back to ./raws (and ./jpegs, ./tiffs) in the
current directory.

Layout produced for a raws folder at  shoot1/raws :
    shoot1/raws/    input ARW files
    shoot1/jpegs/   generated, fed to Claude
    shoot1/tiffs/   generated, get renamed at the end
    shoot1/rename_map.csv

Requires: pip install rawpy imageio exifread anthropic
Set your key:  export ANTHROPIC_API_KEY=...
"""

import sys
import csv
import json
import base64
from pathlib import Path
from collections import Counter

RAW_GLOB = "*.ARW"           # case-sensitive on Linux; see note in stage 0
EXPECTED_COUNTS = (2, 3, 4)  # photos per pair; anything else gets flagged
MODEL = "claude-sonnet-4-5"  # check current docs before a big run


def resolve_paths(raw_arg):
    """Given the raws folder, derive jpegs/, tiffs/ and the CSV beside it."""
    raw_dir = Path(raw_arg).expanduser().resolve() if raw_arg else Path("raws").resolve()
    base = raw_dir.parent
    return {
        "raw": raw_dir,
        "jpeg": base / "jpegs",
        "tiff": base / "tiffs",
        "csv": base / "rename_map.csv",
    }


# ----------------------------------------------------------------------
# Stage 0: convert raws to JPEG + TIFF with a shared stem
# ----------------------------------------------------------------------
def convert(p):
    import rawpy
    import imageio.v3 as iio

    if not p["raw"].is_dir():
        sys.exit(f"Raws folder not found: {p['raw']}")
    p["jpeg"].mkdir(exist_ok=True)
    p["tiff"].mkdir(exist_ok=True)

    raws = sorted(p["raw"].glob(RAW_GLOB))
    if not raws:
        sys.exit(f"No files matching {RAW_GLOB} in {p['raw']}/. "
                 "On Linux the glob is case-sensitive: check .arw vs .ARW.")

    for raw_path in raws:
        with rawpy.imread(str(raw_path)) as raw:
            rgb = raw.postprocess(use_camera_wb=True)  # 8-bit RGB
        stem = raw_path.stem
        iio.imwrite(p["jpeg"] / f"{stem}.jpg", rgb, quality=90)
        iio.imwrite(p["tiff"] / f"{stem}.tif", rgb)
        print(f"converted {stem}")

    print(f"\nDone. {len(raws)} raws -> {p['jpeg']}/ and {p['tiff']}/")


# ----------------------------------------------------------------------
# Stage 1: Claude reads both jar IDs, build rename_map.csv
# ----------------------------------------------------------------------
PROMPT = """This image shows two jars, each labelled with an ID.
Return ONLY a JSON object, no other text:
{"left": "<id on left jar>", "right": "<id on right jar>"}
If an ID is unreadable, use null for that jar."""


def _exif_datetime(path):
    import exifread
    with open(path, "rb") as f:
        tags = exifread.process_file(f, stop_tag="DateTimeOriginal")
    return str(tags.get("EXIF DateTimeOriginal", ""))  # "YYYY:MM:DD HH:MM:SS"


def extract(p):
    import anthropic

    if p["csv"].exists():
        sys.exit(f"{p['csv']} already exists. Move or delete it before re-running "
                 "extract, so a reviewed map is never silently overwritten.")
    if not p["jpeg"].is_dir():
        sys.exit(f"No jpegs folder at {p['jpeg']}. Run convert first.")

    client = anthropic.Anthropic()
    records = []

    for jpg in sorted(p["jpeg"].glob("*.jpg")):
        data = base64.standard_b64encode(jpg.read_bytes()).decode()
        msg = client.messages.create(
            model=MODEL,
            max_tokens=100,
            messages=[{
                "role": "user",
                "content": [
                    {"type": "image",
                     "source": {"type": "base64", "media_type": "image/jpeg", "data": data}},
                    {"type": "text", "text": PROMPT},
                ],
            }],
        )
        try:
            ids = json.loads(msg.content[0].text.strip())
            left, right = ids.get("left"), ids.get("right")
        except (json.JSONDecodeError, AttributeError):
            left, right = None, None

        dt = _exif_datetime(jpg)
        date = dt.split()[0].replace(":", "") if dt else ""
        records.append({
            "old_stem": jpg.stem,
            "left": left, "right": right,
            "date": date,
            "pair_key": f"{left}_{right}_{date}",
            "dt": dt,
        })
        print(f"read {jpg.stem}: left={left} right={right} date={date}")

    # number photos within each pair+date group, ordered by capture time
    records.sort(key=lambda r: (r["pair_key"], r["dt"]))
    counters = {}
    for r in records:
        n = counters.get(r["pair_key"], 0) + 1
        counters[r["pair_key"]] = n
        r["photo"] = n
        r["new_stem"] = f"{r['left']}_{r['right']}_{r['date']}_{n}"

    with open(p["csv"], "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["old_stem", "left", "right", "date", "photo", "new_stem"])
        for r in records:
            w.writerow([r["old_stem"], r["left"], r["right"],
                        r["date"], r["photo"], r["new_stem"]])

    # warn about groups with an unexpected number of photos (likely a misread ID)
    counts = Counter(r["pair_key"] for r in records)
    flagged = {k: c for k, c in counts.items() if c not in EXPECTED_COUNTS}
    print(f"\nWrote {len(records)} rows to {p['csv']}.")
    if flagged:
        print("\nCHECK these groups (unexpected photo count, an ID may be misread):")
        for k, c in flagged.items():
            print(f"  {k}  has {c} photos")
    print(f"\nNow open the CSV and verify it by hand. Then run:\n"
          f"  python renaming_pipeline.py apply {p['raw']} --i-have-reviewed")


# ----------------------------------------------------------------------
# Stage 3: rename TIFFs from the reviewed CSV
# ----------------------------------------------------------------------
def apply(p, reviewed):
    if not reviewed:
        sys.exit("Refusing to run. Open rename_map.csv, check it, then pass "
                 "--i-have-reviewed to confirm.")
    if not p["csv"].exists():
        sys.exit(f"{p['csv']} not found. Run extract first.")

    with open(p["csv"]) as f:
        rows = list(csv.DictReader(f))

    # pre-flight: refuse if any new name is incomplete or would collide
    seen = {}
    problems = []
    for r in rows:
        new = r["new_stem"]
        if "None" in new or new.endswith("_") or new.startswith("None"):
            problems.append(f"  incomplete name from {r['old_stem']}: {new}")
        seen.setdefault(new, []).append(r["old_stem"])
    for new, olds in seen.items():
        if len(olds) > 1:
            problems.append(f"  name collision {new}.tif <- {', '.join(olds)}")
    if problems:
        print("Stopped before renaming anything. Fix these in the CSV first:")
        print("\n".join(problems))
        sys.exit(1)

    renamed = 0
    for r in rows:
        tif = p["tiff"] / f"{r['old_stem']}.tif"
        target = tif.with_stem(r["new_stem"])
        if not tif.exists():
            print(f"skip, no tiff for {r['old_stem']}")
            continue
        if target.exists():
            print(f"skip, target already exists: {target.name}")
            continue
        tif.rename(target)
        renamed += 1
        print(f"{tif.name} -> {target.name}")

    print(f"\nRenamed {renamed} of {len(rows)} TIFFs.")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    cmd = args[0] if args else ""
    raw_arg = args[1] if len(args) > 1 else None
    p = resolve_paths(raw_arg)

    if cmd == "convert":
        convert(p)
    elif cmd == "extract":
        extract(p)
    elif cmd == "apply":
        apply(p, "--i-have-reviewed" in flags)
    else:
        sys.exit(__doc__)