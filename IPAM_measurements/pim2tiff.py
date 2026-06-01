#!/usr/bin/env python3
"""
pim2tiff.py - Convert Walz Imaging-PAM .PIM to a multiframe TIFF, reproducing
ImagingWin's own export, fully headless (no Windows / no ImagingWin).

Validated BIT-EXACT against an ImagingWin TIFF export for the Fv/Fm protocol
(ImagingPAM software 2.41a, MAXI head).

Encoding (reverse-engineered and validated against ground truth):
  * .PIM is a single zlib stream.
  * Decompressed payload = a small TLV text header ("ImagingPAM", version, date,
    time, "Settings: ..."), then 8-bit 640x480 frames, one byte per pixel,
    row-major.
  * ImagingWin multiplies every pixel by 4 on export:  export = stored * 4.
  * Walz frame order: 1=Fo 2=Fm 3=RedAbsorptivity 4=NIRAbsorptivity, then
    F'/Fm' induction pairs (5/6, 7/8, ...). RedAbs/NIRAbs are empty in the
    Fv/Fm protocol.
  * The four base frames are contiguous (1-byte separator between frames). Each
    F'/Fm' pair is preceded by a fixed-size metadata marker (MARKER_LEN). The
    first frame is anchored on the "Settings:" record so the parser adapts to
    different settings-string lengths.

CAVEAT: validated on a 6-frame Fv/Fm file. SCALE=4 and the frame stride are
verified; MARKER_LEN is measured from one file. For longer induction-curve
protocols (more F'/Fm' pairs), spot-check one file against an ImagingWin export
before trusting a whole backlog. Use --marker to override if needed.

Requires: numpy, tifffile
"""
import argparse, sys, zlib
from pathlib import Path
import numpy as np
import tifffile

W, H = 640, 480
FRAME = W * H            # 307200 bytes per 8-bit frame
STRIDE = FRAME + 1       # 1-byte separator between consecutive frames
SCALE = 4                # ImagingWin export multiplier (verified bit-exact)
MARKER_LEN = 729         # metadata block before each F'/Fm' induction pair
DESCRIPTOR = 3           # bytes between the Settings record and the first frame


def _first_frame_offset(payload: bytes) -> int:
    k = payload.find(b"Settings:")
    if k < 0:
        raise ValueError("no 'Settings:' record found; unexpected PIM header")
    tag = payload.rfind(b"\x06", 0, k)      # 0x06 string-record tag
    if tag < 0:
        raise ValueError("malformed Settings record")
    val_len = payload[tag + 1]
    rec_end = tag + 2 + val_len             # first byte after the settings value
    return rec_end + DESCRIPTOR


def parse_pim(path, marker_len: int = MARKER_LEN) -> np.ndarray:
    payload = zlib.decompress(Path(path).read_bytes())
    buf = np.frombuffer(payload, np.uint8)
    n = len(buf)
    pos = _first_frame_offset(payload)

    frames = []
    idx = 0
    while pos + FRAME <= n:
        frames.append(buf[pos:pos + FRAME].reshape(H, W).astype(np.uint16) * SCALE)
        pos += FRAME
        idx += 1
        if pos + FRAME > n:
            break
        # a metadata marker precedes the start of each induction pair (frame
        # index 4, 6, 8, ...); base frames use a 1-byte separator.
        if idx >= 4 and (idx - 4) % 2 == 0:
            pos += marker_len
        else:
            pos += 1
    return np.stack(frames)


def convert(path: Path, outdir: Path, split: bool, marker_len: int):
    frames = parse_pim(path, marker_len)
    outdir.mkdir(parents=True, exist_ok=True)
    multi = outdir / (path.stem + ".tif")
    tifffile.imwrite(multi, frames, photometric="minisblack", metadata={"axes": "ZYX"})
    if split:
        labels = ["Fo", "Fm", "RedAbs", "NIRAbs"]
        for i in range(frames.shape[0]):
            tag = labels[i] if i < 4 else f"frame{i+1:02d}"
            tifffile.imwrite(outdir / f"{path.stem}-{tag}.tif", frames[i],
                             photometric="minisblack")
    nz = [int(np.count_nonzero(frames[i])) for i in range(frames.shape[0])]
    print(f"{path.name}: {frames.shape[0]} frames -> {multi.name} | nonzero/frame={nz}")


def main():
    ap = argparse.ArgumentParser(description="Walz .PIM -> multiframe TIFF (headless).")
    ap.add_argument("inputs", nargs="+", help="one or more .pim files or a folder")
    ap.add_argument("-o", "--outdir", type=Path, required=True)
    ap.add_argument("--split", action="store_true", help="also write per-frame TIFFs")
    ap.add_argument("--marker", type=int, default=MARKER_LEN,
                    help=f"induction-marker length in bytes (default {MARKER_LEN})")
    a = ap.parse_args()
    files = []
    for it in a.inputs:
        p = Path(it)
        files += (sorted(p.glob("*.pim")) + sorted(p.glob("*.PIM"))) if p.is_dir() else [p]
    if not files:
        sys.exit("no .pim files found")
    for f in files:
        try:
            convert(f, a.outdir, a.split, a.marker)
        except Exception as e:
            print(f"FAILED {f.name}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()