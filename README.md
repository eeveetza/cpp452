# C++ Implementation of Recommendation ITU-R P.452

This code repository contains a C++ software implementation of  [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)  with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz.  

## Version
This implementation corresponds to the original reference MATLAB/Octave implementation of [Recommendation ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en) approved by ITU-R Working Party 3M and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx).  

## Folder structure

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`/src/P452.cpp`                | C++ implementation of Recommendation ITU-R P.452-18         |
|`/src/maps/`          | Subfolder containing data maps provided with Recommendation ITU-R P.452 (see [Integrating ITU Digital Products](#integrating-itu-digital-products)). |
|`/tests/`           | Folder containing validation tests against the reference MATLAB/Octave implementation from [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)     |
|`/tests/validate_p452.cpp` | Googletest test that reads every profile/results pair in `/tests/validation_examples/` at run time and checks both the intermediate path-profile parameters (`ae`, `dtot`, `hts`, `hrs`, `theta_t`, `theta_r`, `theta`, `hm`, `hte`, `hre`, `hstd`, `hsrd`, `dlt`, `dlr`, `pathtype`, `dtm`, `dlm`, `b0`, `omega`, `DN`, `N0`) and the final/intermediate transmission-loss terms (`Lbfsg`, `Lb0p`, `Lb0b`, `Ldsph`, `Ld50`, `Ldp`, `Lbs`, `Lba`, `Lb`) against the reference values - mirroring [`validate_p452.m`](https://github.com/eeveetza/p452/blob/main/matlab/validate_p452.m) in the MATLAB/Octave repository. |
|`/tests/validation_examples/` | Copy of the MATLAB repository's `validation_examples/` folder (`profiles/*.csv`, `results/*.csv`). Update this folder whenever the MATLAB repository's copy changes; `validate_p452.cpp` will automatically pick up the new/changed cases the next time the tests run - no regeneration step needed. |
|`/tests/test101.cpp`, `/tests/test102.cpp` | Unit tests for the digital-map bilinear interpolation, independent of the profile-based validation examples above. |


## Integrating ITU Digital Products

This software uses ITU digital products that are integral part of Recommendations. These products must not be reproduced or distributed without explicit written permission from the ITU.

**Download and extract the required maps** to `./src/maps`:
   - From ITU-R P.452-18:
     - `DN50.TXT`
     - `N050.TXT`
- Ensure all files are placed in `./src/maps` before running the code.

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
./build.sh
~~~

The main validation test, `tests/validate_p452.cpp`, reads every profile/results pair under `tests/validation_examples/` at run time (the same files used by [`validate_p452.m`](https://github.com/eeveetza/p452/blob/main/matlab/validate_p452.m) in the MATLAB/Octave repository) and checks both intermediate and final results against them, so keeping `tests/validation_examples/` in sync with the MATLAB repository is enough to keep this test suite current - no code generation or regeneration step is required.

## Software Versions
The code was tested and runs on:
* Apple clang version 4.0.0 (tags/RELEASE_400/final) on x86_64-apple-darwin22.6.0

## Prerequisites
The code uses [googletest](https://github.com/google/googletest) for automatised testing.

## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)

* [MATLAB/Octave Implementation of Recommendation ITU-R P.452](https://github/eeveetza/p452)
