# Helheim ATLAS LiDAR

## Building on Linux

## Create your environment

    conda create -n atlas
    conda activate atlas

## Install PDAL and a compiler

    conda install -c conda-forge pdal compilers ninja cmake eigen -y

## Build

    ./build.sh

## How to run displacement

```bash
./build/atlas <before.las> <after.las> --xshift=<xshift> --yshift=<yshift>
```

- `xshift`/`yshift`: expected bulk displacement in meters (based on Helheim velocity × time between scans)
- `--dumpxy X,Y`: dump debug info for the cell containing UTM coordinate X,Y
- `--dumpij I,J`: dump debug info for grid cell at index I,J

Outputs:
- `<before_stem>.tif`: GeoTIFF (EPSG:32624) with 5 bands
- `vector.svg`: displacement arrow field

| Band | Name | Description |
|------|------|-------------|
| 1 | X | Mean X displacement (m) |
| 2 | Y | Mean Y displacement (m) |
| 3 | Z | Mean Z displacement (m) |
| 4 | BEFORE | Point count before scan |
| 5 | AFTER | Point count after scan |

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
