# C++ Implementation of Recommendation ITU-R P.452

This code repository contains a C++ software implementation of  [Recommendation ITU-R P.452-17](https://www.itu.int/rec/R-REC-P.452/en)  with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz.  



| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`/src/P452.cpp`                | C++ implementation of Recommendation ITU-R P.452-17         |
|`/src/P452.test.cpp`           | Validation tests against the reference MATLAB/Octave implementation from [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)     |



## Function Call

~~~ 
Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_path, Gt, Gr, pol, dct, dcr, DN, N0, press, temp, ha_t, ha_r, dk_t, dk_r);
~~~

## Required input arguments of function `tl_p452`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `f`               | double | GHz   | ~0.1 ≤ `f` ≤ ~50 | Frequency   |
| `p         `      | double  | %     | 0.001 ≤ `p` ≤ 50 | Time percentage for which the calculated basic transmission loss is not exceeded |
| `d`               | `vector<double>` | km    |  0 < `max(d)` ≤ ~10000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | `vector<double>` | m (asl)   |   | Terrain profile heights |
| `zone`           | `vector<int>`   |       | 1 - Coastal Land, 2 - Inland, 3 - Sea             |  Radio-climatic zone types |
| `htg`           | double    | m      |           |  Tx antenna height above ground level |
| `hrg`           | double    | m      |          |  Rx antenna height above ground level |
| `phi_path`           | double    | deg      |   -90 ≤ `phi_path`  ≤ 90          |  Latitude of path center between Tx and Rx stations |
| `Gt`  `Gr`           | double  |   dBi    |           |  Tx/Rx antenna gain in the direction of the horizon saalong the great-circle interference path. |
| `pol`           | int    |       |   `pol`  = 1, 2          |  Polarization of the signal: 1 - horizontal, 2 - vertical |
| `dct`           | double    | km      |   `dct` ≥ 0          |  Distance over land from the Tx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `dcr`           | double    | km      |   `dcr` ≥ 0          |  Distance over land from the Rx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `DN`            | double    | N-units/km      | `DN`> 0           | The average radio-refractivity lapse-rate through the lowest 1 km of the atmosphere at the path-center. It can be derived from an appropriate map (see below).  |
| `N0`           | double    | N-units      |             | The sea-level surface refractivity at the path-center. It can be derived from an appropriate map (see below).|
| `press`           | double    | hPa      |             | Dry air pressure.|
| `temp`           | double    | deg C      |             | Air temperature.|
| `ha_t`           | double    | m      |             | Clutter nominal height at the Tx side |
| `ha_r`           | double    | m      |             | Clutter nominal height at the Rx side |
| `dk_t`           | double    | km      |             | Clutter nominal distance at the Tx side |
| `dk_r`           | double    | km      |             | Clutter nominal distance at the Rx side |

The clutter loss will be computed if `ha_t` > `htg` or `ha_r` > `hrg`.

Note that  `d` needs to be significantly greater than `dk_t` and/or `dk_r` for the clutter model in Recommendation P.452-17 to be applicable.

 
## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lb`    | double | dB    | Basic transmission loss |


## Meteorological Data
The following input arguments related to meteorological data:

* `DN`: The average radio-refractivity lapse-rate through the lowest 1 km of the atmosphere at the path center
* `N0`: The sea-level surface refractivity at the path center

can be derived from the appropriate maps provided with the recommendation, using the following function calls
~~~ 
std::string path = "./R-REC-P/"; // Path to the P452 maps DN50.TXT and N050.TXT
P452DigitalMaps maps = P452DigitalMaps(path); // load P452 maps

double DN50 = maps.GetDN50(lon, lat);  // interpolate the value of DN50 at the point given by (lon, lat)
double N050 = maps.GetN050(lon, lat);  // interpolate the value of N050 at the point given by (lon, lat)
~~~

## Path Center
The path center latitude `phi_path` is either provided with the terrain data, or can be computed using the function `great_circle_path`, which computes the coordinates of the center point (`phime`, `phimn`) 
of the great-circle path between the Tx position (`phite`, `phitn`)  and the Rx position (`phire`, `phirn`) 

~~~
double Re = 6371.0; // average Earth's radius (km)
double dpnt = 0.5;  // midpoint along the great-circle path 
double dgc;         // Tx-Rx distance along great-circle path
double bt2r;        // bearing of the great-circle path between Tx and Rx
great_circle_path(phire, phite, phirn, phitn, Re, dpnt, phime, phimn, bt2r, dgc);
phi_path = phimn;
~~~

## Validation

Validation can be performed by compiling the code and runing the tests as follows

~~~
make
./test
~~~

## Software Versions
The code was tested and runs on:
* Apple clang version 12.0.0 (clang-1200.0.32.29) on x86_64-apple-darwin21.3.0

## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)

* [MATLAB/Octave Implementation of Recommendation ITU-R P.452](https://github/eeveetza/p452)
