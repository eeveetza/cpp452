# C++ Implementation of Recommendation ITU-R P.452

This code repository contains a C++ software implementation of  [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)  with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz.  

## Version
This development branch repository contains a C++ software implementation of  [Recommendation ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en).  

## Folder structure
| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`/src/P452.cpp`                | C++ implementation of Recommendation ITU-R P.452-18         |
|`/src/P452.test.cpp`           | Validation tests against the reference MATLAB/Octave implementation from [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)  (to be published)   |



## Function Call

~~~ 
Lb = tl_p452(maps, f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp);
~~~

## Required input arguments of function `tl_p452`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `maps`           | class `P452DigitalMaps` | |  | Object containing all the digital maps (DN50, N050) necessary for computation |
| `f`               | double | GHz   | ~0.1 ≤ `f` ≤ ~50 | Frequency   |
| `p`               | double  | %     | 0.001 ≤ `p` ≤ 50 | Time percentage for which the calculated basic transmission loss is not exceeded |
| `d`               | `vector<double>` | km    |  0 < `max(d)` ≤ ~10000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | `vector<double>` | m (asl)   |   | Terrain profile heights |
| `g`          | `vector<double>` | m (asl)   |  | Clutter + Terrain profile heights   |
| `zone`           | `vector<int>`   |       | 1 - Coastal Land, 2 - Inland, 3 - Sea             |  Radio-climatic zone types |
| `htg`           | double    | m      |           |  Tx antenna height above ground level |
| `hrg`           | double    | m      |          |  Rx antenna height above ground level |
| `phit_e`           | double    | deg      |   -180 ≤ `phit_e`  ≤ 180            |  Tx Longitude |
| `phit_n`           | double    | deg      |   -90 ≤ `phit_n`  ≤ 90           |  Tx Latitude |
| `phir_e`           | double    | deg      |   -180 ≤ `phir_e`  ≤ 180            |  Rx Longitude |
| `phir_n`           | double    | deg      |   -90 ≤ `phir_n`  ≤ 90           |  Rx Latitude |
| `Gt`  `Gr`           | double  |   dBi    |           |  Tx/Rx antenna gain in the direction of the horizon saalong the great-circle interference path. |
| `pol`           | int    |       |   `pol`  = 1, 2          |  Polarization of the signal: 1 - horizontal, 2 - vertical |
| `dct`           | double    | km      |   `dct` ≥ 0          |  Distance over land from the Tx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `dcr`           | double    | km      |   `dcr` ≥ 0          |  Distance over land from the Rx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `press`           | double    | hPa      |             | Dry air pressure.|
| `temp`           | double    | deg C      |             | Air temperature.|


 
## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lb`    | double | dB    | Basic transmission loss |



## Validation

Validation can be performed by compiling the code and runing the tests as follows

~~~
make
./test
~~~

## Software Versions
The code was tested and runs on:
* Apple clang version 14.0.3 (clang-1403.0.22.14.1) on x86_64-apple-darwin22.5.0


## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)

* [MATLAB/Octave Implementation of Recommendation ITU-R P.452](https://github/eeveetza/p452)
