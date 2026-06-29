# Helheim ATLAS LiDAR

## Building on Linux

## Create your environment

    conda create -n atlas
    conda activate atlas

## Install PDAL and a compiler

    conda install -c conda-forge pdal compilers ninja cmake eigen -y

## Build

```bash
cd displacement
./build.sh
```

`./build.sh` does a full clean build (wipes `build/` and reruns cmake). After
code-only changes you can rebuild incrementally with `cd build && ninja`.

The `--ncc` refinement pass is parallelized with OpenMP. Thread count is
controlled by the `OMP_NUM_THREADS` environment variable (defaults to all
cores); set `OMP_NUM_THREADS=1` for serial execution.

## How to run displacement

```bash
./build/atlas <before.las> <after.las> --xshift=<xshift> --yshift=<yshift> [options]
```

- `xshift`/`yshift`: expected bulk displacement in meters (based on Helheim velocity × time between scans)
- `--topfrac F`: top fraction of points (by height above the tile's fitted plane) used as the surface slice (default 0.5)
- `--gridlen M`: grid cell size in meters for flood-fill shape detection (default 1.0)
- `--minshape N`: minimum shape size in grid cells; smaller shapes dropped before matching (default 1)
- `--slices N`: number of contiguous height bands the top `--topfrac` is split into for shape detection (default 1)
- `--passes N`: grid passes; passes after the first reprocess every cell seeded from its converged neighborhood, roughly N× runtime (default 1)
- `--ncc`: refine each cell's shape-match displacement by normalized cross-correlation of detrended mean-height rasters with parabolic subpixel peak fit; adds bands 14–16 (NCC_X, NCC_Y, NCC_PEAK)
- `--nccradius N`: NCC search radius in raster cells around the shape-match displacement (default 5)
- `--ncclen M`: NCC raster cell size in meters, independent of `--gridlen` (default 2.0; nccradius × ncclen must not exceed the 20 m tile overlap)
- `--tiff DIR`: output directory for GeoTIFF (default: current directory)
- `--geojson DIR`: output directory for shapes GeoJSON (omit to skip)
- `--svg DIR`: output directory for displacement arrow SVG (omit to skip)
- `--dumpxy X,Y`: dump debug info for the cell containing UTM coordinate X,Y
- `--dumpij I,J`: dump debug info for grid cell at index I,J

Outputs:
- `<before_stem>_<after_stem>.tif`: GeoTIFF (EPSG:32624) with 13 bands (16 with `--ncc`)
- `<before_stem>_<after_stem>.geojson`: matched shape polygons (opt-in via `--geojson`)
- `<before_stem>_<after_stem>.svg`: displacement arrow field (opt-in via `--svg`)

> **Note on CRS:** the `atlas` displacement output is **EPSG:32624** (WGS 84 / UTM zone 24N), whereas the published DEM and velocity products in the S3 bucket are **EPSG:3413** (WGS 84 / NSIDC Sea Ice Polar Stereographic North). Reproject before comparing the two.

### Output band layout

| Band | Name | Description |
|------|------|-------------|
| 1 | X | Mean X displacement (m) |
| 2 | Y | Mean Y displacement (m) |
| 3 | Z | Mean Z displacement (m) |
| 4 | BEFORE | Point count before scan |
| 5 | AFTER | Point count after scan |
| 6 | MEDIAN_X | Median X displacement |
| 7 | MEDIAN_Y | Median Y displacement |
| 8 | MEDIAN_Z | Median Z displacement |
| 9 | MATCH_COUNT | Number of matched shape pairs |
| 10 | RMS_RESIDUAL | RMS deviation of matched pairs from mean |
| 11 | MAD_X | Median absolute deviation of X displacement across matched pairs |
| 12 | MAD_Y | Median absolute deviation of Y displacement across matched pairs |
| 13 | MAD_Z | Median absolute deviation of Z displacement across matched pairs |
| 14 | NCC_X | (`--ncc` only) NCC-refined total X displacement (m) |
| 15 | NCC_Y | (`--ncc` only) NCC-refined total Y displacement (m) |
| 16 | NCC_PEAK | (`--ncc` only) peak normalized cross-correlation value |

Bands 9 and 10 are quality indicators — low match count or high RMS residual
signals an unreliable cell. Cells with 1–2 pairs report near-zero RMS/MAD that
doesn't reflect real uncertainty, so filter on MATCH_COUNT (band 9) downstream.

With `--ncc`, bands 1–13 are unchanged (bit-identical to a run without the
flag) — the refinement goes only to bands 14–16 so the two estimates stay
directly comparable per cell. NCC_X/Y are written only where the peak passed
the acceptance gate (peak ≥ 0.5, not on the search-window edge); NCC_PEAK is
written wherever the search produced a scoreable surface.

### Example using ATLAS South COPC files from S3

6-hour scan pair:

```bash
./build/atlas \
    https://atlas-lidar-helheim.s3.us-west-2.amazonaws.com/copc/south/191012_070542.copc.laz \
    https://atlas-lidar-helheim.s3.us-west-2.amazonaws.com/copc/south/201013_130153.copc.laz \
    --xshift=6.4 --yshift=-0.51
```

24-hour scan pair:

```bash
./build/atlas \
    https://atlas-lidar-helheim.s3.us-west-2.amazonaws.com/copc/south/230108_131655.copc.laz \
    https://atlas-lidar-helheim.s3.us-west-2.amazonaws.com/copc/south/230109_131631.copc.laz \
    --xshift=25.6 --yshift=-2.04
```

The xshift/yshift rule of thumb for a 6-hour scan is 6.4 m / -0.51 m; scale linearly with time for other intervals.


## Bucket layout

The dataset is publicly available at `s3://atlas-lidar-helheim` (no credentials required).

```
s3://atlas-lidar-helheim/
├── copc/
│   ├── north/          # ATLAS North Cloud Optimized Point Clouds (.copc.laz)
│   └── south/          # ATLAS South Cloud Optimized Point Clouds (.copc.laz)
├── dem/
│   ├── ATLAS-North/    # ATLAS North digital elevation models (GeoTIFF)
│   └── ATLAS-South/    # ATLAS South digital elevation models (GeoTIFF)
└── velocity/
    ├── ATLAS-North/    # ATLAS North glacier surface velocity rasters (GeoTIFF)
    └── ATLAS-South/    # ATLAS South glacier surface velocity rasters (GeoTIFF)
```

File names follow the convention `YYMMDD_HHMMSS` (UTC). Elevation models are orthometric heights with the 49 m EIGEN-6C4 geoid removed. Velocity GeoTIFFs have 4 bands: absolute horizontal velocity, X velocity, Y velocity, and Z velocity. The DEM and velocity products are in EPSG:3413 (WGS 84 / NSIDC Sea Ice Polar Stereographic North) — note this differs from the EPSG:32624 displacement output of the `atlas` tool above.

To browse the bucket:

```bash
aws s3 ls --no-sign-request s3://atlas-lidar-helheim/
```

Copy a file to your local machine:

```bash
aws s3 cp --no-sign-request s3://atlas-lidar-helheim/dem/ATLAS-North/191215_121054_idw_geoid_rm.cog.tif .
```
