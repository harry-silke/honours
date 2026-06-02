"""
Jar-photo pipeline: ARW -> JPEG + TIFF -> Claude reads jar IDs -> rename TIFFs.

Run as three separate stages. The CSV review between stage 2 and stage 3 is a
HARD STOP: stage 3 refuses to run unless you have looked at the CSV and marked
it approved (see --i-have-reviewed below).

Pass the raws folder as an argument. jpegs/, tiffs/ and rename_map.csv are
created alongside it (in the same parent), so each shoot stays self-contained:

    python renaming_pipeline.py convert /path/to/shoot1/raws
    python renaming_pipeline.py extract /path/to/shoot1/raws
    # ---- check the names, then verify the CSV by hand ----
    python renaming_pipeline.py preview /path/to/shoot1/raws
    python renaming_pipeline.py apply /path/to/shoot1/raws --i-have-reviewed

extract is resumable: re-run it and it skips photos already in the CSV.
preview prints the names apply would create, and renames nothing.

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
JPEG_MAX_DIM = 2400          # longest edge of the JPEG Claude reads (keeps it under 5 MB)
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
    import io
    import rawpy
    import imageio.v3 as iio
    from PIL import Image

    if not p["raw"].is_dir():
        sys.exit(f"Raws folder not found: {p['raw']}")
    p["jpeg"].mkdir(exist_ok=True)
    p["tiff"].mkdir(exist_ok=True)

    raws = [f for f in sorted(p["raw"].glob(RAW_GLOB))
            if not f.name.startswith("._")]   # skip macOS AppleDouble sidecars
    if not raws:
        sys.exit(f"No files matching {RAW_GLOB} in {p['raw']}/. "
                 "On Linux the glob is case-sensitive: check .arw vs .ARW.")

    converted = 0
    for raw_path in raws:
        try:
            with rawpy.imread(str(raw_path)) as raw:
                rgb = raw.postprocess(use_camera_wb=True)  # 8-bit RGB
        except Exception as e:
            # one unreadable file shouldn't kill the whole run
            print(f"SKIP {raw_path.name}: {e}")
            continue
        stem = raw_path.stem

        # full-size TIFF: this is the real data
        iio.imwrite(p["tiff"] / f"{stem}.tif", rgb)

        # downscaled JPEG: only used for Claude to read the labels, stays well
        # under the 5 MB API limit. Copy EXIF so extract can read capture date.
        img = Image.fromarray(rgb)
        img.thumbnail((JPEG_MAX_DIM, JPEG_MAX_DIM))
        exif = b""
        try:
            with open(raw_path, "rb") as f:
                src = Image.open(f)
                exif = src.info.get("exif", b"")
        except Exception:
            pass
        img.save(p["jpeg"] / f"{stem}.jpg", format="JPEG",
                 quality=85, exif=exif)
        converted += 1
        print(f"converted {stem}")

    print(f"\nDone. {converted} of {len(raws)} raws -> {p['jpeg']}/ and {p['tiff']}/")


# ----------------------------------------------------------------------
# Stage 1: Claude reads both jar IDs, build rename_map.csv
# ----------------------------------------------------------------------
PROMPT = """This photo shows petri dishes, each with its clear lid flipped and
resting ABOVE its own dish. A handwritten label in marker is on each lid. The
labels are upside down in the original orientation. To help, a second copy of
the same photo, rotated 180 degrees so the labels read normally, is also provided.

Usually there are TWO lids (left and right, as seen in the original upright-photo
orientation). Sometimes there is only ONE lid.

Each label is a number, sometimes followed by a single letter (the treatment),
for example: 127, 127B, 127C. Read the FULL label including any trailing letter.
Do not drop the letter. Do not add a letter that is not written.

Return ONLY a JSON object, no other text:
{"left": "<full label on left lid>", "right": "<full label on right lid>"}

Rules:
- Two lids present: left label in "left", right label in "right".
- Only ONE lid present: put its label in "left" and use "1000" for "right".
- A label is present but genuinely unreadable: use "1000" for that side.
Never use null. Use the string "1000" for any missing or unreadable side."""


def _exif_datetime(path):
    import exifread
    with open(path, "rb") as f:
        tags = exifread.process_file(f, stop_tag="DateTimeOriginal")
    return str(tags.get("EXIF DateTimeOriginal", ""))  # "YYYY:MM:DD HH:MM:SS"


def _jpeg_views_under_limit(path, max_bytes=5_000_000):
    """Return (original_b64, rotated180_b64) JPEG strings, both downscaled to stay
    under the API limit. The rotated copy makes upside-down labels read normally.
    The file on disk is never modified."""
    import io
    from PIL import Image

    img = Image.open(path)
    if img.mode not in ("RGB", "L"):
        img = img.convert("RGB")

    def encode(im):
        for max_dim in (2400, 1800, 1400, 1000):
            scaled = im.copy()
            scaled.thumbnail((max_dim, max_dim))
            buf = io.BytesIO()
            scaled.save(buf, format="JPEG", quality=85)
            if buf.tell() <= max_bytes:
                return base64.standard_b64encode(buf.getvalue()).decode()
        return base64.standard_b64encode(buf.getvalue()).decode()

    return encode(img), encode(img.rotate(180))


CSV_FIELDS = ["old_stem", "left", "right", "date", "dt"]


def _read_one(client, jpg, raw_dir):
    """Send one photo (plus its rotated copy) to Claude, return a result dict."""
    original_b64, rotated_b64 = _jpeg_views_under_limit(jpg)
    msg = client.messages.create(
        model=MODEL,
        max_tokens=100,
        messages=[{
            "role": "user",
            "content": [
                {"type": "text", "text": "Original photo (labels upside down):"},
                {"type": "image",
                 "source": {"type": "base64", "media_type": "image/jpeg", "data": original_b64}},
                {"type": "text", "text": "Same photo rotated 180 degrees (labels upright):"},
                {"type": "image",
                 "source": {"type": "base64", "media_type": "image/jpeg", "data": rotated_b64}},
                {"type": "text", "text": PROMPT},
            ],
        }],
    )
    reply = msg.content[0].text.strip()
    if reply.startswith("```"):
        reply = reply.strip("`")
        if reply.lstrip().startswith("json"):
            reply = reply.lstrip()[4:]
        reply = reply.strip()
    try:
        ids = json.loads(reply)
        left, right = ids.get("left"), ids.get("right")
    except (json.JSONDecodeError, AttributeError):
        left, right = None, None

    left = "1000" if left in (None, "", "null") else str(left).strip().upper()
    right = "1000" if right in (None, "", "null") else str(right).strip().upper()

    raw_path = raw_dir / f"{jpg.stem}.ARW"
    dt = _exif_datetime(raw_path) if raw_path.exists() else ""
    date = dt.split()[0].replace(":", "") if dt else ""
    return {"old_stem": jpg.stem, "left": left, "right": right,
            "date": date, "dt": dt, "reply": reply}


def extract(p, batch_size=25):
    """Read jar labels into the CSV, one row per photo, written as we go.
    Safe to re-run: photos already in the CSV are skipped, so an interrupted
    run resumes from where it stopped."""
    import anthropic

    if not p["jpeg"].is_dir():
        sys.exit(f"No jpegs folder at {p['jpeg']}. Run convert first.")

    # which photos are already done?
    done = set()
    if p["csv"].exists():
        with open(p["csv"]) as f:
            done = {row["old_stem"] for row in csv.DictReader(f)}

    all_jpgs = [j for j in sorted(p["jpeg"].glob("*.jpg"))
                if not j.name.startswith("._")]   # skip macOS AppleDouble sidecars
    todo = [j for j in all_jpgs if j.stem not in done]
    if not todo:
        print(f"All {len(all_jpgs)} photos already in {p['csv'].name}. Nothing to do.")
        _summarise(p)
        return
    print(f"{len(done)} already done, {len(todo)} to read "
          f"(batches of {batch_size}).")

    client = anthropic.Anthropic()
    new_file = not p["csv"].exists()
    # append mode so a resumed run adds to the existing CSV
    with open(p["csv"], "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if new_file:
            w.writeheader()
        for i in range(0, len(todo), batch_size):
            batch = todo[i:i + batch_size]
            for jpg in batch:
                try:
                    r = _read_one(client, jpg, p["raw"])
                except Exception as e:
                    # one bad photo shouldn't kill the run; leave it for a re-run
                    print(f"ERROR on {jpg.stem}: {e}  (skipped, re-run to retry)")
                    continue
                w.writerow(r)
                f.flush()  # persist immediately so a crash loses nothing
                note = (f"  <- model said: {r['reply'][:80]!r}"
                        if r["left"] == "1000" and r["right"] == "1000" else "")
                print(f"read {r['old_stem']}: left={r['left']} "
                      f"right={r['right']} date={r['date']}{note}")
            print(f"-- batch {i // batch_size + 1} done, CSV saved --")

    _summarise(p)


def _numbered_rows(rows):
    """Given raw per-photo rows, assign photo number within each pair+date group
    (ordered by capture time) and build new_stem. Returns rows with both added."""
    for r in rows:
        r["pair_key"] = f"{r['left']}_{r['right']}_{r['date']}"
    rows.sort(key=lambda r: (r["pair_key"], r["dt"], r["old_stem"]))
    counters = {}
    for r in rows:
        n = counters.get(r["pair_key"], 0) + 1
        counters[r["pair_key"]] = n
        r["photo"] = n
        r["new_stem"] = f"{r['left']}_{r['right']}_{r['date']}_{n}"
    return rows


def _summarise(p):
    """Print the group-count check from the current CSV contents."""
    with open(p["csv"]) as f:
        rows = _numbered_rows(list(csv.DictReader(f)))
    counts = Counter(r["pair_key"] for r in rows)
    flagged = {k: c for k, c in counts.items() if c not in EXPECTED_COUNTS}
    print(f"\nCSV now holds {len(rows)} photos: {p['csv']}")
    if flagged:
        print("\nCHECK these groups (unexpected photo count, a label may be misread):")
        for k, c in flagged.items():
            print(f"  {k}  has {c} photos")
    print(f"\nVerify the CSV, then run:\n"
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
        rows = _numbered_rows(list(csv.DictReader(f)))

    # pre-flight: warn on sentinel values, refuse on collisions
    seen = {}
    sentinels = []
    problems = []
    for r in rows:
        new = r["new_stem"]
        if "1000" in (r["left"], r["right"]):
            sentinels.append(f"  {r['old_stem']} -> {new}  (single/unreadable lid)")
        seen.setdefault(new, []).append(r["old_stem"])
    for new, olds in seen.items():
        if len(olds) > 1:
            problems.append(f"  name collision {new}.tif <- {', '.join(olds)}")
    if sentinels:
        print("Note: these have a 1000 sentinel (filter downstream as planned):")
        print("\n".join(sentinels))
        print()
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
# preview: show the names apply WOULD produce, change nothing
# ----------------------------------------------------------------------
def preview(p):
    if not p["csv"].exists():
        sys.exit(f"{p['csv']} not found. Run extract first.")

    with open(p["csv"]) as f:
        rows = _numbered_rows(list(csv.DictReader(f)))

    seen = {}
    sentinels = 0
    for r in rows:
        print(f"{r['old_stem']}.tif  ->  {r['new_stem']}.tif")
        if "1000" in (r["left"], r["right"]):
            sentinels += 1
        seen.setdefault(r["new_stem"], []).append(r["old_stem"])

    collisions = {n: o for n, o in seen.items() if len(o) > 1}
    counts = Counter(r["pair_key"] for r in rows)
    flagged = {k: c for k, c in counts.items() if c not in EXPECTED_COUNTS}

    print(f"\n{len(rows)} photos. Nothing was renamed (preview only).")
    if sentinels:
        print(f"{sentinels} row(s) carry a 1000 sentinel (single/unreadable lid).")
    if flagged:
        print("\nGroups with an unexpected photo count (a label may be misread):")
        for k, c in flagged.items():
            print(f"  {k}  has {c} photos")
    if collisions:
        print("\nName collisions (apply will refuse until these are fixed):")
        for n, olds in collisions.items():
            print(f"  {n}.tif <- {', '.join(olds)}")


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
    elif cmd == "preview":
        preview(p)
    elif cmd == "apply":
        apply(p, "--i-have-reviewed" in flags)
    else:
        sys.exit(__doc__)