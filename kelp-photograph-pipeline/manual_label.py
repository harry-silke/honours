"""
Manual labelling tool for jar-lid photos.

Pairs with renaming_pipeline.py but doesn't depend on it. Run convert first so
you have a jpegs/ folder, then point this at the same raws path:

    python manual_label.py /path/to/shoot1/raws

The script opens each JPEG in Preview (macOS) and waits for you to type the
labels in the terminal. Each enter saves the row to manual_labels.csv beside
the jpegs and advances to the next photo.

Keys (single character + enter, or just the labels + enter):
    <label> Tab <label> Enter   commit and advance
    Enter on left field          skip this photo for now
    Empty right field            auto-fills "1000" (one lid case)
    b Enter                       go back one photo to re-edit
    q Enter                       quit (progress is already saved)

Resume: re-running picks up where you left off (skips photos already in the CSV).
The output CSV uses the same schema as renaming_pipeline's extract, so you can
run `python renaming_pipeline.py preview <raws>` and `apply` on it directly,
as long as you rename manual_labels.csv to rename_map.csv first (or symlink it).
"""

import sys
import csv
import subprocess
from pathlib import Path

CSV_FIELDS = ["old_stem", "left", "right", "date", "dt"]
CSV_NAME = "manual_labels.csv"


def _exif_datetime(path):
    try:
        import exifread
    except ImportError:
        return ""
    try:
        with open(path, "rb") as f:
            tags = exifread.process_file(f, stop_tag="DateTimeOriginal")
        return str(tags.get("EXIF DateTimeOriginal", ""))
    except Exception:
        return ""


def _open_in_preview(path):
    """Tell Preview to display this image. -g keeps focus in the terminal."""
    subprocess.run(["open", "-g", "-a", "Preview", str(path)], check=False)


def _close_preview():
    """Best-effort close of any open Preview windows when we finish."""
    subprocess.run(
        ["osascript", "-e", 'tell application "Preview" to quit'],
        check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def _normalise(s):
    """Strip and uppercase, but keep '1000' as is."""
    return s.strip().upper() if s else ""


def _load_existing(csv_path):
    """Return a dict of old_stem -> row for everything already labelled, and
    an ordered list of stems in CSV order."""
    if not csv_path.exists():
        return {}, []
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    return {r["old_stem"]: r for r in rows}, [r["old_stem"] for r in rows]


def _append_row(csv_path, row):
    new_file = not csv_path.exists()
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if new_file:
            w.writeheader()
        w.writerow(row)
        f.flush()


def _rewrite_csv(csv_path, ordered_stems, by_stem):
    """Used when going back: rewrite the whole CSV in order, omitting the last."""
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for stem in ordered_stems:
            w.writerow(by_stem[stem])


def _print_recent(history):
    """Print last 5 commits so the user can spot a typo a few rows back."""
    if not history:
        return
    print("\nRecent:")
    for r in history[-5:]:
        print(f"  {r['old_stem']}: {r['left']} / {r['right']}")
    print()


def main(raw_arg):
    raw_dir = Path(raw_arg).expanduser().resolve()
    if not raw_dir.is_dir():
        sys.exit(f"Raws folder not found: {raw_dir}")
    base = raw_dir.parent
    jpeg_dir = base / "jpegs"
    if not jpeg_dir.is_dir():
        sys.exit(f"No jpegs folder at {jpeg_dir}. Run renaming_pipeline.py convert first.")
    csv_path = base / CSV_NAME

    all_jpgs = sorted(jpeg_dir.glob("*.jpg"))
    if not all_jpgs:
        sys.exit(f"No JPEGs in {jpeg_dir}")

    by_stem, order = _load_existing(csv_path)
    done = set(by_stem)
    todo = [j for j in all_jpgs if j.stem not in done]
    history = [by_stem[s] for s in order]

    total = len(all_jpgs)
    print(f"\n{len(done)} already labelled, {len(todo)} to go.")
    print(f"CSV: {csv_path}")
    print(f"Photos open in Preview. Type left<Tab>right<Enter>. "
          f"'b'=back, 'q'=quit, blank right -> 1000.\n")

    i = 0
    while i < len(todo):
        jpg = todo[i]
        _open_in_preview(jpg)
        _print_recent(history)
        print(f"[{len(done) + i + 1}/{total}]  {jpg.name}")

        try:
            raw_input_text = input("  left<Tab>right: ").rstrip("\n")
        except (EOFError, KeyboardInterrupt):
            print("\nQuitting. Progress is saved.")
            break

        cmd = raw_input_text.strip().lower()
        if cmd == "q":
            print("Quitting. Progress is saved.")
            break
        if cmd == "b":
            if not history:
                print("Nothing to go back to.")
                continue
            last = history.pop()
            order.remove(last["old_stem"])
            del by_stem[last["old_stem"]]
            _rewrite_csv(csv_path, order, by_stem)
            todo.insert(i, jpeg_dir / f"{last['old_stem']}.jpg")
            print(f"Stepped back to {last['old_stem']} (was {last['left']} / {last['right']}).")
            continue
        if raw_input_text.strip() == "":
            print("  skipped (no entry recorded)")
            i += 1
            continue

        # parse "left<TAB>right" or "left right" (space fallback)
        if "\t" in raw_input_text:
            parts = raw_input_text.split("\t", 1)
        else:
            parts = raw_input_text.split(None, 1)
        left = _normalise(parts[0])
        right = _normalise(parts[1]) if len(parts) > 1 else ""
        if right == "":
            right = "1000"

        raw_path = raw_dir / f"{jpg.stem}.ARW"
        dt = _exif_datetime(raw_path) if raw_path.exists() else ""
        date = dt.split()[0].replace(":", "") if dt else ""
        row = {"old_stem": jpg.stem, "left": left, "right": right,
               "date": date, "dt": dt}
        _append_row(csv_path, row)
        by_stem[jpg.stem] = row
        order.append(jpg.stem)
        history.append(row)
        i += 1

    _close_preview()
    print(f"\nLabelled {len(history)} of {total} photos. CSV at {csv_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])