"""
GRASS GIS 8.5 batch script: fill nodata in atlas displacement rasters.

Usage:
    python grass_fill_stats.py INPUT*.TIF [options]
    python grass_fill_stats.py /path/to/velocities/*.tif --workers 14 --bands 6

Options:
    --workers N         Parallel GRASS sessions (default: 1)
    --out-dir DIR       Output directory (default: same as input)
    --suffix SUFFIX     Output filename suffix (default: _filled)
    --bands N [N ...]   Bands to fill (default: 6)
    --distance N        Search radius in cells (default: 3)
    --mode MODE         wmean|mean|median|mode (default: wmean)
    --power F           IDW power (default: 2.0)
    --cells N           Min data cells in radius (default: 8)
    --keep              Preserve existing values, only fill NULLs
    --skip-existing     Skip output files that already exist (safe to resume)
    --min-match       Null out cells with MATCH_COUNT < N before filling

"""

import json
import os
import re
import shutil
import subprocess
import sys
import argparse
import logging
import tempfile
import threading
from datetime import datetime
from pathlib import Path

from tqdm import tqdm

# Single-token line the inner GRASS session prints to stdout after each file;
# the outer process reads these to advance one shared tqdm progress bar.
PROGRESS_TOKEN = "__PROGRESS__"

def _find_grass_env_bin() -> str:
    """Locate the bin/ dir holding grass, gdal_translate and gdalinfo.

    Order: $GRASS_ENV_BIN override, then the dir of `grass` on PATH (i.e. the
    active conda env), so the same script runs on any machine with grass_env
    activated. Returns "" if not found; main() reports the error.
    """
    override = os.environ.get("GRASS_ENV_BIN")
    if override:
        return override
    grass = shutil.which("grass")
    return os.path.dirname(grass) if grass else ""


GRASS_ENV_BIN = _find_grass_env_bin()
GRASS_BIN = os.path.join(GRASS_ENV_BIN, "grass")
GDAL_TRANSLATE = os.path.join(GRASS_ENV_BIN, "gdal_translate")
GDALINFO = os.path.join(GRASS_ENV_BIN, "gdalinfo")


def sanitize(name: str) -> str:
    """Convert a file stem to a valid GRASS map name."""
    return re.sub(r"[^A-Za-z0-9]", "_", name)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("inputs", nargs="+", metavar="INPUT.TIF")
    p.add_argument("--workers", type=int, default=1,
                   help="Parallel GRASS sessions (default: 1)")
    p.add_argument("--out-dir", default=None,
                   help="Write filled TIFs here (default: same dir as input)")
    p.add_argument("--suffix", default="_filled",
                   help="Appended to the input stem for the output filename")
    p.add_argument("--bands", nargs="+", type=int, default=[6], metavar="N",
                   help="Bands to fill (default: 6)")
    p.add_argument("--distance", type=int, default=3,
                   help="r.fill.stats search radius in cells (default: 3)")
    p.add_argument("--mode", choices=["wmean", "mean", "median", "mode"],
                   default="wmean")
    p.add_argument("--power", type=float, default=2.0,
                   help="IDW power coefficient (default: 2.0)")
    p.add_argument("--cells", type=int, default=8,
                   help="Minimum data cells within search radius (default: 8)")
    p.add_argument("--min-match", type=int, default=None, metavar="N",
                   help="Null out cells with MATCH_COUNT < N before filling")
    p.add_argument("--max-rms", type=float, default=None, metavar="F",
                   help="Null out cells with RMS_RESIDUAL > F before filling")
    p.add_argument("--keep", dest="no_keep", action="store_false",
                   help="Preserve existing values, only fill NULLs")
    p.add_argument("--skip-existing", action="store_true",
                   help="Skip files whose output already exists")
    p.add_argument("--_inner", action="store_true", help=argparse.SUPPRESS)
    return p.parse_args(argv) 


# ---------------------------------------------------------------------------
# Per-file processing (runs inside a GRASS session)
# ---------------------------------------------------------------------------

def output_path(tif: Path, args) -> Path:
    dest_dir = Path(args.out_dir) if args.out_dir else tif.parent
    return dest_dir / f"{tif.stem}{args.suffix}.tif"


def process_tif(gs, tif: Path, args):
    out = output_path(tif, args)

    if args.skip_existing and out.exists():
        logging.info("[skip] %s", tif.name)
        return

    stem = sanitize(tif.stem)
    fill_bands = sorted(args.bands)
    keep_flag = "" if args.no_keep else "k"
    filled_tifs = {}

    src_info = json.loads(subprocess.run(
        [GDALINFO, "-json", str(tif)],
        capture_output=True, text=True, check=True,
    ).stdout)
    src_bands = {b["band"]: b for b in src_info.get("bands", [])}
    src_descriptions = {
        n: b.get("description", f"BAND_{n}") for n, b in src_bands.items()
    }
    src_nodata = {
        n: b.get("noDataValue") for n, b in src_bands.items()
    }

    with tempfile.TemporaryDirectory(prefix="grass_fill_") as tmpdir:
        use_quality = args.min_match is not None or args.max_rms is not None
        map_match = map_rms = None

        if use_quality:
            for qband, qname in ((9, "match"), (10, "rms")):
                qfile = os.path.join(tmpdir, f"q{qband}.tif")
                qnodata = src_nodata.get(qband)
                qcmd = [GDAL_TRANSLATE, "-b", str(qband), str(tif), qfile]
                if qnodata is not None:
                    qcmd[2:2] = ["-a_nodata", str(qnodata)]
                subprocess.run(qcmd, capture_output=True, check=True)
                qmap = f"{stem}_q{qname}"
                gs.run_command("r.in.gdal", input=qfile, output=qmap,
                               flags="ek", overwrite=True)
                if qname == "match":
                    map_match = qmap
                else:
                    map_rms = qmap

        for b in fill_bands:
            label = src_descriptions.get(b, str(b))
            logging.info("[%s] band %s (%s) ...", tif.name, b, label)

            nodata = src_nodata.get(b)

            # Extract just this band as a temp single-band TIF
            extracted = os.path.join(tmpdir, f"b{b}_in.tif")
            extract_cmd = [GDAL_TRANSLATE, "-b", str(b), str(tif), extracted]
            if nodata is not None:
                extract_cmd[2:2] = ["-a_nodata", str(nodata)]
            subprocess.run(extract_cmd, capture_output=True, check=True)

            map_in = f"{stem}_b{b}_in"
            map_out = f"{stem}_b{b}_out"
            filled = os.path.join(tmpdir, f"b{b}_filled.tif")

            gs.run_command("r.in.gdal", input=extracted, output=map_in,
                           flags="ek", overwrite=True)
            gs.run_command("g.region", raster=map_in)

            if use_quality:
                map_masked = f"{stem}_b{b}_masked"
                conditions = []
                if args.min_match is not None:
                    conditions.append(f"{map_match} >= {args.min_match}")
                if args.max_rms is not None:
                    conditions.append(f"{map_rms} <= {args.max_rms}")
                condition = " && ".join(conditions)
                gs.run_command("r.mapcalc",
                               expression=f"{map_masked} = if({condition}, {map_in}, null())",
                               overwrite=True)
                fill_input = map_masked
            else:
                fill_input = map_in

            gs.run_command(
                "r.fill.stats",
                input=fill_input, output=map_out,
                distance=args.distance, mode=args.mode,
                power=args.power, cells=args.cells,
                flags=keep_flag, overwrite=True,
            )
            out_gdal_kwargs = dict(
                input=map_out, output=filled,
                format="GTiff", type="Float32", overwrite=True,
                flags="c",  # suppress color table metadata (avoids harmless GTiff ERROR 6)
            )
            if nodata is not None:
                out_gdal_kwargs["nodata"] = nodata
            gs.run_command("r.out.gdal", **out_gdal_kwargs)
            remove_maps = f"{map_in},{map_out}"
            if use_quality:
                remove_maps += f",{map_masked}"
            gs.run_command("g.remove", type="raster",
                           name=remove_maps, flags="f", quiet=True)

            filled_tifs[b] = filled

        if use_quality and map_match:
            gs.run_command("g.remove", type="raster",
                           name=f"{map_match},{map_rms}", flags="f", quiet=True)

        out.parent.mkdir(parents=True, exist_ok=True)

        width, height = src_info["size"]
        gt_str = ", ".join(str(v) for v in src_info["geoTransform"])
        srs_wkt = src_info["coordinateSystem"]["wkt"]

        # Build VRT with Description tags, one output band per filled band
        lines = [
            f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}">',
            f'  <SRS>{srs_wkt}</SRS>',
            f'  <GeoTransform>{gt_str}</GeoTransform>',
        ]
        for out_band, src_band in enumerate(fill_bands, start=1):
            desc = src_descriptions.get(src_band, f"BAND_{src_band}")
            nodata = src_nodata.get(src_band)
            nodata_line = f"    <NoDataValue>{nodata}</NoDataValue>" if nodata is not None else ""
            lines += [
                f'  <VRTRasterBand dataType="Float32" band="{out_band}">',
                f'    <Description>{desc}</Description>',
                *([nodata_line] if nodata_line else []),
                f'    <SimpleSource>',
                f'      <SourceFilename relativeToVRT="0">{filled_tifs[src_band]}</SourceFilename>',
                f'      <SourceBand>1</SourceBand>',
                f'    </SimpleSource>',
                f'  </VRTRasterBand>',
            ]
        lines.append("</VRTDataset>")

        vrt = os.path.join(tmpdir, "combined.vrt")
        Path(vrt).write_text("\n".join(lines))
        subprocess.run(
            [GDAL_TRANSLATE, "-of", "GTiff", "-co", "COMPRESS=LZW",
             vrt, str(out)],
            capture_output=True, check=True,
        )

    logging.info("wrote %s", out)


# ---------------------------------------------------------------------------
# Worker launch helpers
# ---------------------------------------------------------------------------

def make_inner_cmd(script_path: str, batch_files: list, args) -> list:
    cmd = [GRASS_BIN, "--tmp-project", "EPSG:32624", "--exec",
           "python3", script_path, "--_inner"]
    cmd += batch_files
    cmd += ["--bands"] + [str(b) for b in args.bands]
    cmd += ["--distance", str(args.distance), "--mode", args.mode]
    cmd += ["--power", str(args.power), "--cells", str(args.cells)]
    if args.min_match is not None:
        cmd += ["--min-match", str(args.min_match)]
    if args.max_rms is not None:
        cmd += ["--max-rms", str(args.max_rms)]
    if not args.no_keep:
        cmd += ["--keep"]
    if args.skip_existing:
        cmd += ["--skip-existing"]
    if args.out_dir:
        cmd += ["--out-dir", str(args.out_dir)]
    cmd += ["--suffix", args.suffix]
    return cmd


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def write_run_header(logf, args, files, workers, script_path):
    """Record the full run configuration at the top of the log so every
    output directory can be traced back to the exact settings that made it."""
    out_dir = Path(args.out_dir).resolve() if args.out_dir else "same as input"
    input_dirs = sorted({str(Path(f).parent) for f in files})
    lines = [
        "=" * 72,
        f"grass_fill_stats run  {datetime.now().isoformat(timespec='seconds')}",
        f"host:          {os.uname().nodename}",
        f"cwd:           {os.getcwd()}",
        f"script:        {script_path}",
        f"python:        {sys.executable}",
        f"grass:         {GRASS_BIN}",
        f"input files:   {len(files)}",
        f"input dirs:    {', '.join(input_dirs)}",
        f"out dir:       {out_dir}",
        f"suffix:        {args.suffix}",
        f"workers:       {workers} (requested {args.workers})",
        "-- r.fill.stats settings --",
        f"bands:         {' '.join(str(b) for b in sorted(args.bands))}",
        f"distance:      {args.distance} cells",
        f"mode:          {args.mode}",
        f"power:         {args.power}",
        f"cells:         {args.cells}",
        f"keep (-k):     {not args.no_keep}  "
        f"({'only NULLs filled' if not args.no_keep else 'all cells smoothed'})",
        f"min-match:     {args.min_match}",
        f"max-rms:       {args.max_rms}",
        f"skip-existing: {args.skip_existing}",
        "-- command (inputs abbreviated) --",
        " ".join(
            [Path(sys.argv[0]).name, f"<{len(files)} input files>"]
            + [a for a in sys.argv[1:] if a not in set(args.inputs)]
        ),
        "=" * 72,
    ]
    logf.write("\n".join(lines) + "\n")
    # Flush before workers start: they inherit this file as stderr and must
    # append after the header, not interleave with Python's write buffer.
    logf.flush()


def _pump_progress(proc, pbar, lock):
    """Read a worker's stdout; advance the shared bar on each progress token."""
    for line in proc.stdout:
        if line.strip() == PROGRESS_TOKEN:
            with lock:
                pbar.update(1)
    proc.stdout.close()


def main():
    args = parse_args()

    missing = [b for b in (GRASS_BIN, GDAL_TRANSLATE, GDALINFO)
               if not os.access(b, os.X_OK)]
    if missing:
        sys.exit(
            "GRASS/GDAL executables not found: " + ", ".join(missing) + "\n"
            "Activate the GRASS conda env (conda activate grass_env) or set "
            "GRASS_ENV_BIN to the env's bin/ directory."
        )

    if not args._inner and "GISBASE" not in os.environ:
        files = []
        for f in args.inputs:
            p = Path(f).resolve()
            if p.is_file():
                files.append(str(p))
            else:
                print(f"WARNING: skipping missing file: {f}", file=sys.stderr)

        if not files:
            sys.exit("No input files found.")

        workers = min(args.workers, len(files))
        script_path = os.path.abspath(__file__)

        # Quiet GRASS module chatter (progress %, info) in the workers; warnings
        # and errors still come through. Workers' stderr is captured to a log
        # file so the terminal shows only the tqdm bar.
        os.environ["GRASS_VERBOSE"] = "0"
        log_dir = Path("./log")
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"grass_fill_stats_{datetime.now():%Y%m%d_%H%M%S}.log"

        # Round-robin split so each worker gets an interleaved mix of files
        # (avoids one worker getting all the slow files if sorted by size).
        batches = [files[i::workers] for i in range(workers) if files[i::workers]]

        lock = threading.Lock()
        procs, threads = [], []
        with open(log_path, "w") as logf, tqdm(
            total=len(files), desc="r.fill.stats", unit="file"
        ) as pbar:
            write_run_header(logf, args, files, workers, script_path)
            for batch in batches:
                p = subprocess.Popen(
                    make_inner_cmd(script_path, batch, args),
                    stdout=subprocess.PIPE, stderr=logf, text=True,
                )
                t = threading.Thread(
                    target=_pump_progress, args=(p, pbar, lock), daemon=True
                )
                t.start()
                procs.append(p)
                threads.append(t)

            for t in threads:
                t.join()
            codes = [p.wait() for p in procs]

        failed = sum(1 for c in codes if c != 0)
        if failed:
            print(
                f"WARNING: {failed} worker(s) exited with errors. See {log_path}",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"All workers finished. Log: {log_path}")
        return

    # --- Inner mode: running inside a GRASS session ---
    # Diagnostics go to stderr (captured to the worker log file by the parent);
    # stdout is reserved for progress tokens that drive the parent's tqdm bar.
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(process)d] %(levelname)s %(message)s",
        stream=sys.stderr,
        force=True,
    )

    # `grass --exec` always sets GISBASE; fall back to the env's lib/grass*
    # (version-agnostic) only if this was launched some other way.
    gisbase = os.environ.get("GISBASE")
    if not gisbase:
        candidates = sorted(
            Path(GRASS_ENV_BIN).parent.joinpath("lib").glob("grass[0-9]*")
        )
        if not candidates:
            sys.exit("GISBASE not set and no lib/grass* found next to "
                     f"{GRASS_ENV_BIN}")
        gisbase = str(candidates[-1])
    sys.path.insert(0, os.path.join(gisbase, "etc", "python"))
    import grass.script as gs  # noqa: PLC0415

    for raw in args.inputs:
        tif = Path(raw).resolve()
        if not tif.is_file():
            logging.warning("skipping missing file: %s", tif)
        else:
            try:
                process_tif(gs, tif, args)
            except subprocess.CalledProcessError as e:
                stderr = e.stderr.decode() if e.stderr else ""
                logging.error("[%s]: %s", tif.name, stderr)
            except Exception as e:
                logging.error("[%s]: %s", tif.name, e)
        # Signal completion of this file (success, skip, or error) to the parent.
        print(PROGRESS_TOKEN, flush=True)

    logging.info("Done.")


if __name__ == "__main__":
    main()
