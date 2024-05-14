#include <sstream>
#include "P452.hpp"

// trigonometric functions with arguments in degrees

double sind(double theta_deg) {

    double theta_rad = theta_deg * PI / 180.0;
    double y = sin(theta_rad);
    return y;
}

double cosd(double theta_deg) {

    double theta_rad = theta_deg * PI / 180.0;
    double y = cos(theta_rad);
    return y;
}

double atan2d(double y, double x) {

    double res = atan2(y, x) * 180.0 / PI;
    return res;

}

double asind(double y) {

    double x = asin(y) * 180.0 / PI;
    return x;

}

void great_circle_path(double Phire, double Phite, double Phirn, double Phitn, double Re, double dpnt, double & Phipnte, double & Phipntn, double & Bt2r, double & dgc) {
    //great_circle_path Great-circle path calculations according to P.2001-4 Attachment H
    //   This function computes the great-circle intermediate points on the
    //   radio path as defined in ITU-R P.2001-4 Attachment H
    //
    //     Input parameters:
    //     Phire   -   Receiver longitude, positive to east (deg)
    //     Phite   -   Transmitter longitude, positive to east (deg)
    //     Phirn   -   Receiver latitude, positive to north (deg)
    //     Phitn   -   Transmitter latitude, positive to north (deg)
    //     Re      -   Average Earth radius (km)
    //     dpnt    -   Distance from the transmitter to the intermediate point (km)
    //
    //     Output parameters:
    //     Phipnte -   Longitude of the intermediate point (deg)
    //     Phipntn -   Latitude of the intermediate point (deg)
    //     Bt2r    -   Bearing of the great-circle path from Tx towards the Rx (deg)
    //     dgc     -   Great-circle path length (km)
    //
    //
    //     Rev   Date        Author                          Description
    //     -------------------------------------------------------------------------------
    //     v0   03MAY22      Ivica Stevanovic OFCOM          initial C++ version

    //// H.2 Path length and bearing

    // Difference (deg) in longitude between the terminals (H.2.1)

    double Dlon = Phire - Phite;

    // Calculate quantity r (H.2.2)

    double r = sind(Phitn) * sind(Phirn) + cosd(Phitn) * cosd(Phirn) * cosd(Dlon);

    // Calculate the path length as the angle subtended at the center of
    // average-radius Earth (H.2.3)

    double Phid = acos(r);  // radians

    // Calculate the great-circle path length (H.2.4)

    dgc = Phid * Re;  // km

    // Calculate the quantity x1 (H.2.5a)

    double x1 = sind(Phirn) - r * sind(Phitn);

    // Calculate the quantity y1 (H.2.5b)

    double y1 = cosd(Phitn) * cosd(Phirn) * sind(Dlon);

    // Calculate the bearing of the great-circle path for Tx to Rx (H.2.6)

    Bt2r = Phire;

    if (abs(x1) < 1e-9 && abs(y1) < 1e-9) {
        Bt2r = Phire;
    } else {
        Bt2r = atan2d(y1, x1);
    }

    //// H.3 Calculation of intermediate path point

    // Calculate the distance to the point as the angle subtended at the center
    // of average-radius Earth (H.3.1)

    double Phipnt = dpnt  / Re;  //radians

    // Calculate quantity s (H.3.2)

    double s = sind(Phitn) * cos(Phipnt) + cosd(Phitn) * sin(Phipnt) * cosd(Bt2r);

    // The latitude of the intermediate point is now given by (H.3.3)

    Phipntn = asind(s); // degs

    // Calculate the quantity x2 (H.3.4a)

    double x2 = cos(Phipnt) - s * sind(Phitn);

    // Calculate the quantity y2 (H.3.4b)

    double y2 = cosd(Phitn) * sin(Phipnt) * sind(Bt2r);

    // Calculate the longitude of the intermediate point Phipnte (H.3.5)

    Phipnte = Bt2r;

    if (x2 < 1e-9 && y2 < 1e-9) {
        Phipnte = Bt2r;
    } else {
        Phipnte = Phite + atan2d(y2, x2);
    }

    
}


void earth_rad_eff(double DN, double & ae, double & ab) {
      //earth_rad_eff Median value of the effective Earth radius
      //     [ae, ab] = earth_rad_eff(DN)
      //     This function computes the median value of the effective earth
      //     radius, and the effective Earth radius exceeded for beta0% of time
      //     as defined in ITU-R P.452-17.
      //
      //     Input arguments:
      //     DN      -   the average radiorefractive index lapse-rate through the
      //                 lowest 1 km of the atmosphere (N-units/km)
      //
      //     Output arguments:
      //     ae      -   the median effective Earth radius (km)
      //     ab      -   the effective Earth radius exceeded for beta0 % of time
      //
      //     Example:   
      //     [ae, ab] = earth_rad_eff(DN)
      //
      //     Rev   Date        Author                          Description
      //     -------------------------------------------------------------------------------
      //     v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++


      double k50 = 157 / (157 - DN);     // (5)

      ae = 6371 * k50;          // (6a)

      double kbeta = 3;

      ab = 6371 * kbeta;        // (6b)

}


double inv_cum_norm( double x ) {
      //inv_cum_norm approximation to the inverse cummulative normal distribution
      //   I = inv_cum_norm( x )
      //   This function implements an approximation to the inverse cummulative
      //   normal distribution function for x <= 0.5 as defined in Attachment 3 to
      //   Annex 1 of the ITU-R P.452-17
      //
      //     Rev   Date        Author                          Description
      //     -------------------------------------------------------------------------------
      //     v0    31MAR22     Ivica Stevanovic, OFCOM         Initial version in C++

      if (x < 0.000001) {
      x = 0.000001;
      }

      //if (x > 0.5) {
      //    warning('This function is defined for arguments not larger than 0.5');
      //}

      double tx = sqrt(-2 * log(x));    // eq (172a)

      double C0 = 2.515516698;        // eq (172c)
      double C1 = 0.802853;           // eq (172d)
      double C2 = 0.010328;           // eq (172e)
      double D1 = 1.432788;           // eq (172f)
      double D2 = 0.189269;           // eq (172g)
      double D3 = 0.001308;           // eq (172h)

      double ksi = ((C2 * tx + C1) * tx + C0) / (((D3 * tx + D2) * tx + D1) * tx + 1);  // eq (172b)

      return ksi - tx;            // eq (172)

}


void smooth_earth_heights(std::vector<double> &  d, std::vector<double> & hi, double htg, double hrg, double ae, double f, 
                                          double & hst, double & hsr, double & hstd, double & hsrd, double & hte, double & hre, double & hm, double & dlt, double & dlr, double & theta_t, double & theta_r, double & theta_tot, int & pathtype) {
      
      //smooth_earth_heights smooth-Earth effective antenna heights according to ITU-R P.452-17
      // [hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f)
      // This function derives smooth-Earth effective antenna heights according to
      // Sections 4 and 5 of the Annex 2 of ITU-R P.452-17
      //
      // Input parameters:
      // d          -   vector of terrain profile distances from Tx [0,rDist] (km)
      // hi         -   vector of terrain profile heigths amsl (m)
      // htg, hrg  -   Tx and Rx antenna heights above ground level (m)
      // ae         -   median effective Earth's radius (c.f. Eq (6a))
      // f       -   frequency (GHz)
      //
      // Output parameters:
      //
      // [0] hst, [1] hsr     -   Tx and Rx antenna heigts of the smooth-Earth surface amsl (m)
      // [2] hstd, [3] hsrd   -   Tx and Rx effective antenna heigts for the diffraction model (m)
      // [4] hte, [5] hre     -   Tx and Rx terminal effective heights for the ducting/layer reflection model (m)
      // [6] hm               -   The terrain roughness parameter (m)
      // [7]     dlt          -   interfering antenna horizon distance (km)
      // [8]     dlr          -   Interfered-with antenna horizon distance (km)
      // [9]     theta_t      -   Interfering antenna horizon elevation angle (mrad)
      // [10]    theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
      // [11]    theta_tot    -   Angular distance (mrad)
      // [12]    pathtype     -   1 = 'los', 2 = 'transhorizon'
      //
      // Rev   Date        Author                          Description
      // -------------------------------------------------------------------------------
      // v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++


      int n = d.size();

      double dtot = d[n - 1];

      //Tx and Rx antenna heights above mean sea level amsl (m)
      double hts = hi[0] + htg;
      double hrs = hi[n - 1] + hrg;

      // Section 5.1.6.2


      double v1 = 0;
      for (int ii = 1; ii < n; ii++) {
            v1 = v1 + (d[ii] - d[ii - 1]) * (hi[ii] + hi[ii - 1]);  // Eq (161)
      }
      double v2 = 0;
      for (int ii = 1; ii < n; ii++) {
            v2 = v2 + (d[ii] - d[ii - 1]) * (hi[ii] * (2 * d[ii] + d[ii - 1]) + hi[ii - 1] * (d[ii] + 2 * d[ii - 1]));  // Eq (162)
      }

      hst = (2 * v1 * dtot - v2) / (dtot * dtot);       // Eq (163)
      hsr = (v2 - v1 * dtot) / (dtot * dtot);          // Eq (164)

      // Section 5.1.6.3

      double hobs = hi[1] - (hts * (dtot - d[1]) + hrs * d[1]) / dtot;  // Eq (165d)
      double alpha_obt = hobs / d[1];
      double alpha_obr = hobs / (dtot - d[1]);

      for (int ii = 2; ii < n - 1; ii++) {

            double hh = hi[ii] - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot;  // Eq (165d)

            hobs = std::max(hh, hobs);                 // Eq (165a)

            alpha_obt = std::max(alpha_obt, hh / d[ii]); // Eq (165b)

            alpha_obr = std::max(alpha_obr, hh / (dtot - d[ii])); // Eq (165c)
      }

      // Calculate provisional values for the Tx and Rx smooth surface heights

      double gt = alpha_obt / (alpha_obt + alpha_obr);         // Eq (166e)
      double gr = alpha_obr / (alpha_obt + alpha_obr);         // Eq (166f)

      double hstp = 0;
      double hsrp = 0;

      if (hobs <= 0) {
            hstp = hst;                                 // Eq (166a)
            hsrp = hsr;                                 // Eq (166b)
      } else {
            hstp = hst - hobs * gt;                       // Eq (166c)
            hsrp = hsr - hobs * gr;                       // Eq (166d)
      }

      // calculate the final values as required by the diffraction model

      hstd = 0;
      hsrd = 0;

      if (hstp >= hi[0]) {
            hstd = hi[0];                                // Eq (167a)
      } else {
            hstd = hstp;                                // Eq (167b)
      }

      if (hsrp > hi[n - 1]) {
            hsrd = hi[n - 1];                              // Eq (167c)
      } else {
            hsrd = hsrp;                                // Eq (167d)
      }

      // Interfering antenna horizon elevation angle and distance


      theta_t = 1000 * atan((hi[1] - hts) / (1000 * d[1]) - d[1] / (2 * ae));  // Eq (152)
      int lt = 0;

      double theta;

      for (int ii = 2; ii < n - 1; ii++) {

            theta = 1000 * atan((hi[ii] - hts) / (1000 * d[ii]) - d[ii] / (2 * ae));  // Eq (152)

            if (theta > theta_t) {   // Eq (154)
                  theta_t = theta;
                  lt = ii;
            }


      }

      double theta_td = 1000 * atan((hrs - hts) / (1000 * dtot) - dtot / (2 * ae));  // Eq (153)
      double theta_rd = 1000 * atan((hts - hrs) / (1000 * dtot) - dtot / (2 * ae));  // Eq (156a)

      pathtype = 1; //los

      if (theta_t > theta_td) {   // Eq (150): test for the trans-horizon path
            pathtype = 2; //transhorizon
      }

      dlt = d[lt];                             // Eq (155)

      // Interfered-with antenna horizon elevation angle and distance

      theta_r = 1000 * atan((hi[1] - hrs) / (1000 * (dtot - d[1])) - (dtot - d[1]) / (2 * ae));  // Eq (157)

      int lr = 0;

      for (int ii = 2; ii < n - 1; ii++) {

            theta = 1000 * atan((hi[ii] - hrs) / (1000 * (dtot - d[ii])) - (dtot - d[ii]) / (2 * ae));  // Eq (157)

            if (theta > theta_r) {   // Eq (158)
                  theta_r = theta;
                  lr = ii;
            }


      }

      dlr = dtot - d[lr];                            // Eq (158)

      if (pathtype == 1) {

      theta_t = theta_td;
      theta_r = theta_rd;

      double lambda = 0.2998 / f;
      double Ce = 1 / ae;


      lt = 0;
      double numax = (hi[1] + 500 * Ce * d[1] * (dtot - d[1]) - (hts * (dtot - d[1]) + hrs * d[1]) / dtot) * sqrt(0.002 * dtot / (lambda * d[1] * (dtot - d[1])));

      for (int ii = 2; ii < n - 1; ii++) {


            double nu = (hi[ii] + 500 * Ce * d[ii] * (dtot - d[ii]) - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot) * sqrt(0.002 * dtot / (lambda * d[ii] * (dtot - d[ii])));

            if (nu > numax) {
                  numax = nu;
                  lt = ii;
            }
      }


      dlt = d[lt];
      dlr = dtot - dlt;

      lr = 0;

      for (int ii = 2; ii < n - 1; ii++) {

            if (dlr < dtot - d[ii]) {
                  lr = ii;
            }

      }
      }

      // Angular distance

      theta_tot = 1e3 * dtot / ae + theta_t + theta_r;         // Eq (159)


      // Section 5.1.6.4 Ducting/layer-reflection model

      // Calculate the smooth-Earth heights at transmitter and receiver as
      // required for the roughness factor

      hst = std::min(hst, hi[0]);                           // Eq (168a)
      hsr = std::min(hsr, hi[n - 1]);                       // Eq (168b)

      // Slope of the smooth-Earth surface

      double m = (hsr - hst) / dtot;                          // Eq (169)

      // The terminal effective heigts for the ducting/layer-reflection model

      hte = htg + hi[0] - hst;                       // Eq (170)
      hre = hrg + hi[n - 1] - hsr;

      hm = hi[lt] - (hst + m * d[lt]);

      for (int ii = lt + 1; ii <= lr; ii++) {

      hm = std::max(hm, hi[ii] - (hst + m * d[ii]));
      }


}


 double longest_cont_dist(std::vector<double> & d, std::vector <int> & zone, int zone_r) {
        //longest_cont_dist Longest continuous path belonging to the zone_r
        //     dm = longest_cont_dist(d, zone, zone_r)
        //     This function computes the longest continuous section of the
        //     great-circle path (km) for a given zone_r
        //
        //     Input arguments:
        //     d       -   vector of distances in the path profile
        //     zone    -   vector of zones in the path profile
        //     zone_r  -   reference zone for which the longest continuous section
        //                 is computed
        //
        //     Output arguments:
        //     dm      -   the longest continuous section of the great-circle path (km) for a given zone_r
        //
        //     Example:
        //     dm = longest_cont_dist(d, zone, zone_r)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++

        double dm = 0;
        double dmc = 0;
        double delta;
        int n = d.size();
        if (zone_r == 12) {
            for (int i = 0; i < n; i++)
            {
                if ((zone[i] == 1) || (zone[i] == 2))
                {
                    if (i == 0) {
                        delta = (d[1] - d[0]) / 2.0;
                    } else if (i == n - 1) {
                        delta = (d[n - 1] - d[n - 2]) / 2.0;
                    } else {
                        delta = (d[i + 1] - d[i - 1]) / 2.0;
                    }

                    dmc = dmc + delta;
                    dm = std::max(dm, dmc);
                } else {
                    dmc = 0;
                }
            }
            return dm;
        }
        else
        {
            for (int i = 0; i < n; i++) {
                if (zone[i] == zone_r)
                {
                    if (i == 0) {
                        delta = (d[1] - d[0]) / 2.0;
                    } else if (i == n - 1) {
                        delta = (d[n - 1] - d[n - 2]) / 2.0;
                    } else {
                        delta = (d[i + 1] - d[i - 1]) / 2.0;
                    }

                    dmc = dmc + delta;
                    dm = std::max(dm, dmc);
                }
                else
                {
                    dmc = 0;
                }
            }
            return dm;
        }
    }

 double path_fraction_sea(vector<double> & d, vector<int> & zone, int zone_r) {
        //path_fraction Path fraction belonging to a given zone_r
        //     omega = path_fraction(d, zone, zone_r)
        //     This function computes the path fraction belonging to a given zone_r
        //     of the great-circle path (km)
        //
        //     Input arguments:
        //     d       -   vector of distances in the path profile
        //     zone    -   vector of zones in the path profile
        //     zone_r  -   reference zone for which the fraction is computed
        //
        //     Output arguments:
        //     omega   -   path fraction belonging to the given zone_r
        //
        //     Example:
        //     omega = path_fraction(d, zone, zone_r)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++

        double dm = 0;

        int n = d.size();
        double delta;

        for (int i = 0; i < n; i++) {
            if (zone[i] == zone_r) {
                if (i == 0) {
                    delta = (d[1] - d[0]) / 2.0;
                } else if (i == n - 1) {
                    delta = (d[n - 1] - d[n - 2]) / 2.0;
                } else {
                    delta = (d[i + 1] - d[i - 1]) / 2.0;
                }

                dm = dm + delta;
            }
        }


        return dm / (d[n - 1] - d[0]);
    }


double beta0(double phi, double dtm, double dlm) {
        ////
        //     This function computes the time percentage for which refractive index
        //     lapse-rates exceeding 100 N-units/km can be expected in the first 100
        //     m of the lower atmosphere
        //     as defined in ITU-R P.452-17.
        //
        //     Input arguments:
        //     phi     -   path centre latitude (deg)
        //     dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
        //     dlm     -   the longest continuous inland section of the great-circle path (km)
        //
        //     Output arguments:
        //     b0      -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
        //
        //     Example:
        //     b0 = beta0(phi, dtm, dlm)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    21MAR22     Ivica Stevanovic, OFCOM         First implementation in C++


        double tau = 1 - exp(-(4.12 * 1e-4 * pow(dlm, 2.41)));       // (3a)

        double mu1 = pow(pow(10, -dtm / (16 - 6.6 * tau)) + pow(10, -5 * (0.496 + 0.354 * tau)), 0.2);

        mu1 = std::min(mu1, 1.0);
        double b0 = 0;
        if (abs(phi) <= 70) {
            double mu4 = pow(10, ((-0.935 + 0.0176 * abs(phi)) * log10(mu1)) );   // (4)
            b0 = pow(10, (-0.015 * abs(phi) + 1.67)) * mu1 * mu4;           // (2)
        } else {
            double mu4 = pow(10, (0.3 * log10(mu1)));                            // (4)
            b0 = 4.17 * mu1 * mu4;                                    // (2)
        }

        return b0;
}


    double  p676_ga_ver11(double f, double p, double rho, double T) {
        //p676_ga_ver11 Specific attenuation due to dry air and water vapour
        // p676_ga_ver11(f, p, rho, T)
        // This function computes the specific attenuation due to dry air and water vapour,
        // at frequencies up to 1 000 GHz for different values of of pressure, temperature
        // and humidity by means of a summation of the individual resonance lines from
        // oxygen and water vapour according to ITU-R P.676-11
        //
        // Input parameters:
        // f       -   Frequency (GHz)
        // p       -   Dry air pressure (hPa)
        // rho     -   Water vapor density (g/m^3)
        // T       -   Temperature (K)
        //
        // Output parameters:
        // g  -   sum of specific attenuations due to dry air and water vapour
        //
        //
        // Rev   Date        Author                          Description
        // ----------------------------------------------------------------------------------------------
        // v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++
        


        //// spectroscopic data for oxigen
        //       f0        a1    a2     a3   a4     a5     a6
        double  oxigen[44][7] = {
                {50.474214, 0.975, 9.651, 6.690, 0.0, 2.566, 6.850},
                {50.987745, 2.529, 8.653, 7.170, 0.0, 2.246, 6.800},
                {51.503360, 6.193, 7.709, 7.640, 0.0, 1.947, 6.729},
                {52.021429, 14.320, 6.819, 8.110, 0.0, 1.667, 6.640},
                {52.542418, 31.240, 5.983, 8.580, 0.0, 1.388, 6.526},
                {53.066934, 64.290, 5.201, 9.060, 0.0, 1.349, 6.206},
                {53.595775, 124.600, 4.474, 9.550, 0.0, 2.227, 5.085},
                {54.130025, 227.300, 3.800, 9.960, 0.0, 3.170, 3.750},
                {54.671180, 389.700, 3.182, 10.370, 0.0, 3.558, 2.654},
                {55.221384, 627.100, 2.618, 10.890, 0.0, 2.560, 2.952},
                {55.783815, 945.300, 2.109, 11.340, 0.0, -1.172, 6.135},
                {56.264774, 543.400, 0.014, 17.030, 0.0, 3.525, -0.978},
                {56.363399, 1331.800, 1.654, 11.890, 0.0, -2.378, 6.547},
                {56.968211, 1746.600, 1.255, 12.230, 0.0, -3.545, 6.451},
                {57.612486, 2120.100, 0.910, 12.620, 0.0, -5.416, 6.056},
                {58.323877, 2363.700, 0.621, 12.950, 0.0, -1.932, 0.436},
                {58.446588, 1442.100, 0.083, 14.910, 0.0, 6.768, -1.273},
                {59.164204, 2379.900, 0.387, 13.530, 0.0, -6.561, 2.309},
                {59.590983, 2090.700, 0.207, 14.080, 0.0, 6.957, -0.776},
                {60.306056, 2103.400, 0.207, 14.150, 0.0, -6.395, 0.699},
                {60.434778, 2438.000, 0.386, 13.390, 0.0, 6.342, -2.825},
                {61.150562, 2479.500, 0.621, 12.920, 0.0, 1.014, -0.584},
                {61.800158, 2275.900, 0.910, 12.630, 0.0, 5.014, -6.619},
                {62.411220, 1915.400, 1.255, 12.170, 0.0, 3.029, -6.759},
                {62.486253, 1503.000, 0.083, 15.130, 0.0, -4.499, 0.844},
                {62.997984, 1490.200, 1.654, 11.740, 0.0, 1.856, -6.675},
                {63.568526, 1078.000, 2.108, 11.340, 0.0, 0.658, -6.139},
                {64.127775, 728.700, 2.617, 10.880, 0.0, -3.036, -2.895},
                {64.678910, 461.300, 3.181, 10.380, 0.0, -3.968, -2.590},
                {65.224078, 274.000, 3.800, 9.960, 0.0, -3.528, -3.680},
                {65.764779, 153.000, 4.473, 9.550, 0.0, -2.548, -5.002},
                {66.302096, 80.400, 5.200, 9.060, 0.0, -1.660, -6.091},
                {66.836834, 39.800, 5.982, 8.580, 0.0, -1.680, -6.393},
                {67.369601, 18.560, 6.818, 8.110, 0.0, -1.956, -6.475},
                {67.900868, 8.172, 7.708, 7.640, 0.0, -2.216, -6.545},
                {68.431006, 3.397, 8.652, 7.170, 0.0, -2.492, -6.600},
                {68.960312, 1.334, 9.650, 6.690, 0.0, -2.773, -6.650},
                {118.750334, 940.300, 0.010, 16.640, 0.0, -0.439, 0.079},
                {368.498246, 67.400, 0.048, 16.400, 0.0, 0.000, 0.000},
                {424.763020, 637.700, 0.044, 16.400, 0.0, 0.000, 0.000},
                {487.249273, 237.400, 0.049, 16.000, 0.0, 0.000, 0.000},
                {715.392902, 98.100, 0.145, 16.000, 0.0, 0.000, 0.000},
                {773.839490, 572.300, 0.141, 16.200, 0.0, 0.000, 0.000},
                {834.145546, 183.100, 0.145, 14.700, 0.0, 0.000, 0.000}
        };


        //// spectroscopic data for water-vapor
        //            f0       b1    b2    b3   b4   b5   b6
        double vapor[35][7] = {
                {22.235080, .1079, 2.144, 26.38, .76, 5.087, 1.00},
                {67.803960, .0011, 8.732, 28.58, .69, 4.930, .82},
                {119.995940, .0007, 8.353, 29.48, .70, 4.780, .79},
                {183.310087, 2.273, .668,  29.06, .77, 5.022, .85},
                {321.225630, .0470, 6.179, 24.04, .67, 4.398, .54},
                {325.152888, 1.514, 1.541, 28.23, .64, 4.893, .74},
                {336.227764, .0010, 9.825, 26.93, .69, 4.740, .61},
                {380.197353, 11.67, 1.048, 28.11, .54, 5.063, .89},
                {390.134508, .0045, 7.347, 21.52, .63, 4.810, .55},
                {437.346667, .0632, 5.048, 18.45, .60, 4.230, .48},
                {439.150807, .9098, 3.595, 20.07, .63, 4.483, .52},
                {443.018343, .1920, 5.048, 15.55, .60, 5.083, .50},
                {448.001085, 10.41, 1.405, 25.64, .66, 5.028, .67},
                {470.888999, .3254, 3.597, 21.34, .66, 4.506, .65},
                {474.689092, 1.260, 2.379, 23.20, .65, 4.804, .64},
                {488.490108, .2529, 2.852, 25.86, .69, 5.201, .72},
                {503.568532, .0372, 6.731, 16.12, .61, 3.980, .43},
                {504.482692, .0124, 6.731, 16.12, .61, 4.010, .45},
                {547.676440, .9785, 0.158, 26.00, .70, 4.500, 1.00},
                {552.020960, .1840, 0.158, 26.00, .70, 4.500, 1.00},
                {556.935985, 497.0, 0.159, 30.86, .69, 4.552, 1.00},
                {620.700807, 5.015, 2.391, 24.38, .71, 4.856, .68},
                {645.766085, .0067, 8.633, 18.00, .60, 4.000, .50},
                {658.005280, .2732, 7.816, 32.10, .69, 4.140, 1.00},
                {752.033113, 243.4, 0.396, 30.86, .68, 4.352, .84},
                {841.051732, .0134, 8.177, 15.90, .33, 5.760, .45},
                {859.965698, .1325, 8.055, 30.60, .68, 4.090, .84},
                {899.303175, .0547, 7.914, 29.85, .68, 4.530, .90},
                {902.611085, .0386, 8.429, 28.65, .70, 5.100, .95},
                {906.205957, .1836, 5.110, 24.08, .70, 4.700, .53},
                {916.171582, 8.400, 1.441, 26.73, .70, 5.150, .78},
                {923.112692, .0079, 10.293, 29.0, .70, 5.000, .80},
                {970.315022, 9.009, 1.919, 25.50, .64, 4.940, .67},
                {987.926764, 134.6, 0.257, 29.85, .68, 4.550, .90},
                {1780.000000, 17506., .952, 196.3, 2.00, 24.15, 5.00}
        };


        double theta = 300.0 / T;

        double e = rho * T / 216.7;        // equation (4)
        double g_0 = 0;
        double SumSiFi = 0.0;
        double d = 5.6e-4 * (p + e) * pow(theta, 0.8);                            // equation (9)
        double Ndf = f * p * theta * theta * (6.14e-5 / (d * (1 + pow(f / d, 2.0))) + 1.4e-12 * p * pow(theta, 1.5) / (1 + 1.9e-5 * pow(f, 1.5)));       // equation (8)

        //// Oxigen computation
        for (int i = 0; i < 44; i++) {
            double fi = oxigen[i][0];

            double Si = oxigen[i][1] * 1e-7 * p * pow(theta, 3) * exp(oxigen[i][2] * (1.0 - theta));       // equation (3)

            double df = oxigen[i][3] * 1e-4 * (p * pow(theta, (0.8 - oxigen[i][4])) + 1.1 * e * theta);  // equation (6a)

            // Doppler broadening

            df = sqrt(df * df + 2.25e-6);                                   // equation (6b)

            double delta = (oxigen[i][5] + oxigen[i][6] * theta) * 1e-4 * (p + e) * pow(theta, 0.8);     // equation (7)

            double Fi = f / fi * ((df - delta * (fi - f)) / (pow(fi - f, 2.0) + df * df) + (df - delta * (fi + f)) / (pow(fi + f, 2.0) + df * df));        // equation (5)


            // specific attenuation due to dry air (oxygen, pressure induced nitrogen
            // and non-resonant Debye attenuation), equations (1-2)


            SumSiFi = SumSiFi + Si * Fi;

        }


        g_0 = 0.182 * f * (SumSiFi + Ndf);

        //// vapor computation


        SumSiFi = 0.0;

        for (int i = 0; i < 35; i++) {
            double fi = vapor[i][0];

            double Si = vapor[i][1] * 1e-1 * e * pow(theta, 3.5) * exp(vapor[i][2] * (1.0 - theta));      // equation (3)

            double df = vapor[i][3] * 1e-4 * (p * pow(theta, vapor[i][4]) + vapor[i][5] * e * pow(theta, vapor[i][6]));     // equation (6a)

            // doppler broadening

            df = 0.535 * df + sqrt(0.217 * df * df + 2.1316e-12 * fi * fi / theta); // equation (6b)

            double delta = 0;                                                           // equation (7)

            double Fi = f / fi * ((df - delta * (fi - f)) / (pow(fi - f, 2.0) + df * df) + (df - delta * (fi + f)) / (pow(fi + f, 2.0) + df * df));              // equation (5)

            SumSiFi = SumSiFi + Si * Fi;

        }
        // specific attenuation due to water vapour, equations (1-2)

        double g_w = 0.182 * f * SumSiFi;

        return g_0 + g_w;

    }


    void pl_los(double d, double f, double p, double b0, double w, double temp, double press, double dlt, double dlr, 
                double & Lbfsg, double & Lb0p, double & Lb0b) {
        //pl_los Line-of-sight transmission loss according to ITU-R P.452-17
        //     This function computes line-of-sight transmission loss (including short-term effects)
        //     as defined in ITU-R P.452-17.
        //
        //     Input parameters:
        //     d       -   Great-circle path distance (km)
        //     f       -   Frequency (GHz)
        //     p       -   Required time percentage(s) for which the calculated basic
        //                 transmission loss is not exceeded (//)
        //     b0      -   Point incidence of anomalous propagation for the path
        //                 central location (//)
        //     w       -   Fraction of the total path over water (//)
        //     temp    -   Temperature (degrees C)
        //     press   -   Dry air pressure (hPa)
        //     dlt     -   For a transhorizon path, distance from the transmit antenna to
        //                 its horizon (km). For a LoS path, each is set to the distance
        //                 from the terminal to the profile point identified as the Bullington
        //                 point in the diffraction method for 50// time
        //     dlr     -   For a transhorizon path, distance from the receive antenna to
        //                 its horizon (km). The same note as for dlt applies here.
        //
        //     Output parameters:
        //     Lbfsg   -   Basic transmission loss due to free-space propagation and
        //                 attenuation by atmospheric gases
        //     Lb0p    -   Basic transmission loss not exceeded for time percentage, p//, due to LoS propagation
        //     Lb0b    -   Basic transmission loss not exceedd for time percentage, b0//, due to LoS propagation
        //
        //     Example:
        //     [Lbfsg, Lb0p, Lb0b] = pl_los(d, f, p, b0, w, temp, press, dlt, dlr)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++

        double T = temp + 273.15;

        // water vapor density
        double rho = 7.5 + 2.5 * w;  // (9a)

        // compute specific attenuation due to dry air and water vapor:

        double g = p676_ga_ver11(f, press, rho, T);

        double Ag = g * d;  //(9)

        // Basic transmission loss due to free-space propagation and attenuation
        // by atmospheric gases
        Lbfsg = 92.4 + 20.0 * log10(f) + 20.0 * log10(d) + Ag;  // (8)

            // Corrections for multipath and focusing effects at p and b0
        double Esp = 2.6 * (1 - exp(-0.1 * (dlt + dlr))) * log10(p / 50);   //(10a)
        double Esb = 2.6 * (1 - exp(-0.1 * (dlt + dlr))) * log10(b0 / 50);  //(10b)

        // Basic transmission loss not exceeded for time percentage p// due to
        // LoS propagation
        Lb0p = Lbfsg + Esp;    //(11)

        // Basic transmission loss not exceeded for time percentage b0// due to
        // LoS propagation
        Lb0b = Lbfsg + Esb;    //(12)


    }


double tl_anomalous(double dtot, double dlt, double dlr, double dct, double dcr, double dlm, double hts, double hrs, double hte, double hre, double hm, double theta_t, double theta_r, double f, double p, double temp, double press, double omega, double ae, double b0) {
        //tl_anomalous Basic transmission loss due to anomalous propagation according to P.452-17
        //   Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, ae, b0)
        //
        //   This function computes the basic transmission loss occuring during
        //   periods of anomalous propagation (ducting and layer reflection)
        //   as defined in ITU-R P.452-17 (Section 4.4)
        //
        //     Input parameters:
        //     dtot         -   Great-circle path distance (km)
        //     dlt          -   interfering antenna horizon distance (km)
        //     dlr          -   Interfered-with antenna horizon distance (km)
        //     dct, dcr     -   Distance over land from the transmit and receive
        //                      antennas tothe coast along the great-circle interference path (km).
        //                      Set to zero for a terminal on a ship or sea platform
        //     dlm          -   the longest continuous inland section of the great-circle path (km)
        //     hts, hrs     -   Tx and Rx antenna heights aobe mean sea level amsl (m)
        //     hte, hre     -   Tx and Rx terminal effective heights for the ducting/layer reflection model (m)
        //     hm           -   The terrain roughness parameter (m)
        //     theta_t      -   Interfering antenna horizon elevation angle (mrad)
        //     theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
        //     f            -   frequency expressed in GHz
        //     p            -   percentage of time
        //     temp         -   Temperature (deg C)
        //     press        -   Dry air pressure (hPa)
        //     omega        -   fraction of the total path over water
        //     ae           -   the median effective Earth radius (km)
        //     b0           -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
        //
        //     Output parameters:
        //     Lba    -   the basic transmission loss due to anomalous propagation
        //               (ducting and layer reflection)
        //
        //     Example:
        //     Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, b0)
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    17NOV16     Ivica Stevanovic, OFCOM         Initial version in Java


        // empirical correction to account for the increasing attenuation with
        // wavelength inducted propagation (47a)

        double Alf = 0;

        if (f < 0.5) {
            Alf = 45.375 - 137.0 * f + 92.5 * f * f;
        }

        // site-shielding diffraction losses for the interfering and interfered-with
        // stations (48)

        double theta_t1 = theta_t - 0.1 * dlt;    // eq (48a)
        double theta_r1 = theta_r - 0.1 * dlr;

        double Ast = 0;
        double Asr = 0;

        if (theta_t1 > 0) {
            Ast = 20 * log10(1 + 0.361 * theta_t1 * sqrt(f * dlt)) + 0.264 * theta_t1 * pow(f, 1.0 / 3.0);
        }

        if (theta_r1 > 0) {
            Asr = 20 * log10(1 + 0.361 * theta_r1 * sqrt(f * dlr)) + 0.264 * theta_r1 * pow(f, 1.0 / 3.0);
        }

        // over-sea surface duct coupling correction for the interfering and
        // interfered-with stations (49) and (49a)

        double Act = 0;
        double Acr = 0;

        if (dct <= 5) {
            if (dct <= dlt) {
                if (omega >= 0.75) {
                    Act = -3 * exp(-0.25 * dct * dct) * (1 + tanh(0.07 * (50 - hts)));
                }
            }
        }

        if (dcr <= 5) {
            if (dcr <= dlr) {
                if (omega >= 0.75) {
                    Acr = -3 * exp(-0.25 * dcr * dcr) * (1 + tanh(0.07 * (50 - hrs)));
                }
            }
        }

        // specific attenuation (51)

        double gamma_d = 5e-5 * ae * pow(f, 1.0 / 3.0);

        // angular distance (corrected where appropriate) (52-52a)

        theta_t1 = theta_t;
        theta_r1 = theta_r;

        if (theta_t > 0.1 * dlt) {
            theta_t1 = 0.1 * dlt;
        }

        if (theta_r > 0.1 * dlr) {
            theta_r1 = 0.1 * dlr;
        }

        double theta1 = 1e3 * dtot / ae + theta_t1 + theta_r1;

        double dI = min(dtot - dlt - dlr, 40.0);   // eq (56a)

        double mu3 = 1;

        if (hm > 10) {

            mu3 = exp(-4.6e-5 * (hm - 10) * (43 + 6 * dI));  // eq (56)

        }

        double tau = 1 - exp(-(4.12e-4 *pow(dlm, 2.41)));       // eq (3a)

        double epsilon = 3.5;

        double alpha = -0.6 - epsilon * 1e-9 * pow(dtot, 3.1) * tau;   // eq (55a)

        if (alpha < -3.4) {
             alpha = -3.4;
        }

        // correction for path geometry:

        double mu2 = pow(500 / ae * dtot * dtot / pow(sqrt(hte) + sqrt(hre), 2.0), alpha);

        if (mu2 > 1.0){
            mu2 = 1.0;
        }
        double beta = b0 * mu2 * mu3;      // eq (54)

        double Gamma = 1.076 / pow(2.0058 - log10(beta), 1.012) * exp(-(9.51 - 4.8 * log10(beta) + 0.198 * pow(log10(beta), 2.0)) * 1e-6 * pow(dtot, 1.13));

        // time percentage variablity (cumulative distribution):

        double Ap = -12 + (1.2 + 3.7e-3 * dtot) * log10(p / beta) + 12 * pow(p / beta, Gamma);  // eq (53)

        // time percentage and angular-distance dependent losses within the
        // anomalous propagation mechanism

        double Adp = gamma_d * theta1 + Ap;   // eq (50)

        // gaseous absorbtion derived from equation (9) using rho = 3 g/m^3 for the
        // whole path length

        // water vapor density

        double rho = 7.5 + 2.5 * omega;

        double T = temp + 273.15;

        // compute specific attenuation due to dry air and water vapor:
        double g = p676_ga_ver11(f, press, rho, T);

        double Ag = g * dtot;  //(9)

        // total of fixed coupling losses (except for local clutter losses) between
        // the antennas and the anomalous propagation structure within the
        // atmosphere (47)

        double Af = 102.45 + 20 * log10(f) + 20 * log10(dlt + dlr) + Alf + Ast + Asr + Act + Acr;

        // total basic transmission loss occuring during periods of anomalaous
        // propagation

        return Af + Adp + Ag;

    }


  double tl_tropo(double dtot, double theta, double f, double p, double temp, double press, double N0, double Gt, double Gr ) {
        //tl_tropo Basic transmission loss due to troposcatterer to P.452-17
        //   Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr )
        //
        //   This function computes the basic transmission loss due to troposcatterer
        //   not exceeded for p% of time
        //   as defined in ITU-R P.452-17 (Section 4.3)
        //
        //     Input parameters:
        //     dtot    -   Great-circle path distance (km)
        //     theta   -   Path angular distance (mrad)
        //     f       -   frequency expressed in GHz
        //     p       -   percentage of time
        //     temp    -   Temperature (deg C)
        //     press   -   Dry air pressure (hPa)
        //     N0      -   path centre sea-level surface refractivity derived from Fig. 6
        //     Gt,Gr   -   Antenna gain in the direction of the horizon along the
        //                 great-circle interference path (dBi)
        //
        //     Output parameters:
        //     Lbs    -   the basic transmission loss due to troposcatterer
        //                not exceeded for p% of time
        //
        //     Example:
        //     Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr )
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         Initial version


        //// Body of function

        double T = temp + 273.15;

        // Frequency dependent loss

        double Lf = 25 * log10(f) - 2.5 * pow(log10(f / 2.0), 2.0);    // eq (45a)

        // aperture to medium coupling loss (dB)

        double Lc = 0.051 * exp(0.055 * (Gt + Gr));             // eq (45b)

        // gaseous absorbtion derived from equation (9) using rho = 3 g/m^3 for the
        // whole path length

        double rho = 3;

        // compute specific attenuation due to dry air and water vapor:
        double g = p676_ga_ver11(f, press, rho, T);

        double Ag = g * dtot;  //(9)

        // the basic transmission loss due to troposcatter not exceeded for any time
        // percentage p, below 50// is given

        return 190 + Lf + 20 * log10(dtot) + 0.573 * theta - 0.15 * N0 + Lc + Ag - 10.1 * pow(-log10(p / 50.0), 0.7);
    }

double knife_edge_loss (double x){

      double y = 6.9 + 20.0 * log10(sqrt(pow(x - 0.1, 2.0) + 1.0) + x - 0.1);

      return y;

}


double dl_se_ft_inner(double repsr, double rsigma, double rDist, double rhte, double rhre, double radft, double rFreq, int rpol) {
        //   This function computes the first-term part of Spherical-Earth diffraction
        //   loss exceeded for p% time for antenna heights
        //   as defined in Sec. 4.2.2.1 of the ITU-R P.452-17, equations (30-37)

        //     Input parameters:
        //     repsr    -   Relative permittivity
        //     rsigma   -   Conductivity (S/m)
        //     rDist    -   Great-circle path distance (km)
        //     rhte     -   Effective height of interfering antenna (m)
        //     rhre     -   Effective height of interfered-with antenna (m)
        //     radft    -   effective Earth radius (km)
        //     rf       -   frequency (GHz)
        //     rpol     -   polarization (1-horizontal, 2-vertical)
        //
        //     Output parameters:
        //     rLdft   -   The first-term spherical-Earth diffraction loss not exceeded for p% time
        //                implementing equations (30-37)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++


        double rLdft = 0.0;

        // Normalized factor for surface admittance for horizontal (1) and vertical
        // (2) polarizations

        double rK = 0.036 * pow((radft * rFreq), (-1.0 / 3.0)) * pow(pow((repsr - 1), 2.0) + pow((18 * rsigma / rFreq), 2.0), (-1.0 / 4.0));   // Eq (30a)

        if (rpol == 2.0) {// vertical polarization

            rK = rK * sqrt(pow(repsr, 2.0) + pow(18 * rsigma / rFreq, 2.0));       // Eq (30b)

        }

        // Earth ground/polarization parameter

        double rbeta_dft = (1 + 1.6 * pow(rK, 2.0) + 0.67 * pow(rK, 4.0)) / (1 + 4.5 * pow(rK, 2.0) + 1.53 * pow(rK, 4.0));  // Eq (31)

        // Normalized distance

        double rX = 21.88 * rbeta_dft * pow(rFreq / (radft * radft), (1.0 / 3.0)) * rDist;          // Eq (32)

        // Normalized transmitter and receiver heights

        double rYt = 0.9575 * rbeta_dft * pow((rFreq * rFreq / radft), (1.0 / 3.0)) * rhte;       // Eq (33a)

        double rYr = 0.9575 * rbeta_dft * pow((rFreq * rFreq / radft), (1.0 / 3.0)) * rhre;       // Eq (33b)

        // Calculate the distance term given by:

        double rFx = 0.0;

        if (rX >= 1.6) {
            rFx = 11 + 10 * log10(rX) - 17.6 * rX;
        } else {
            rFx = -20 * log10(rX) - 5.6488 * pow(rX, 1.425);     // Eq (34)
        }

        double rBt = rbeta_dft * rYt;             // Eq (36b)

        double rBr = rbeta_dft * rYr;              // Eq (36b)

        double rGYt = 0.0;
        double rGYr = 0.0;

        if (rBt > 2) {
            rGYt = 17.6 * pow(rBt - 1.1, 0.5) - 5 * log10(rBt - 1.1) - 8;
        } else {
            rGYt = 20 * log10(rBt + 0.1 * pow(rBt, 3.0));
        }

        if (rBr > 2) {
            rGYr = 17.6 * pow(rBr - 1.1, 0.5) - 5 * log10(rBr - 1.1) - 8;
        } else {
            rGYr = 20 * log10(rBr + 0.1 * pow(rBr, 3.0));
        }

        double rLimit = 2 + 20 * log10(rK);

        rGYt = max(rGYt, rLimit);
        rGYr = max(rGYr, rLimit);

        rLdft = -rFx - rGYt - rGYr;

        return rLdft;
    }


double dl_se_ft(double d, double hte, double hre, double adft, double f, double omega, int pol) {
        //dl_se_ft First-term part of spherical-Earth diffraction according to ITU-R P.452-17
        //   This function computes the first-term part of Spherical-Earth diffraction
        //   loss exceeded for p// time for antenna heights
        //   as defined in Sec. 4.2.2.1 of the ITU-R P.452-17
        //
        //     Input parameters:
        //     d       -   Great-circle path distance (km)
        //     hte     -   Effective height of interfering antenna (m)
        //     hre     -   Effective height of interfered-with antenna (m)
        //     adft    -   effective Earth radius (km)
        //     f       -   frequency (GHz)
        //     omega   -   fraction of the path over sea
        //     pol     -   polarization (1-horizontal, 2- vertical)
        //
        //     Output parameters:
        //     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p// time
        //                Ldft(1) is for the horizontal polarization
        //                Ldft(2) is for the vertical polarization
        //
        //     Example:
        //     Ldft = dl_se_ft(d, hte, hre, adft, f, omega)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++


        // First-term part of the spherical-Earth diffraction loss over land

        double epsr;
        double sigma;


        //// Body of function

        // First-term part of the spherical-Earth diffraction loss over land

        epsr = 22;
        sigma = 0.003;

        double Ldft_land = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f, pol);

        // First-term part of the spherical-Earth diffraction loss over sea

        epsr = 80;
        sigma = 5;

        double Ldft_sea = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f, pol);


        // First-term spherical diffraction loss

        return omega * Ldft_sea + (1 - omega) * Ldft_land;      // Eq (29)

    }

double dl_se(double rDist, double rhte, double rhre, double rap, double rFreq, double omega, double rpol) {
        //   This function computes the Spherical-Earth diffraction loss exceeded
        //   for p% time for antenna heights rhte and rhre (m)
        //   as defined in Sec. 4.2.2 of the ITU-R P.452-17
        //   Only land paths are implemented
        //   Input parameters:
        //     rDist    -   Great-circle path distance (km)
        //     rhte     -   Effective height of interfering antenna (m)
        //     rhre     -   Effective height of interfered-with antenna (m)
        //     rap      -   the effective Earth radius in kilometers
        //     rFreq    -   frequency expressed in GHz
        //     omega   -   the fraction of the path over sea
        //     rpol     -   polarization 1- horizontal, 2- vertical
        //
        //     Output parameters:
        //     rLdsph   -   The spherical-Earth diffraction loss not exceeded for p% time
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         Initial version in C++

        double rLambda = 0.2998 / rFreq;
        double rLdsph = 0.0;

        // Calculate the marginal LoS distance for a smooth path

        double rdlos = sqrt(2 * rap) * (sqrt(0.001 * rhte) + sqrt(0.001 * rhre));    // Eq (23)

        if (rDist >= rdlos) {
            // calculate diffraction loss Ldft using the method in Sec. 4.2.2.1 for
            // adft = ap and set Ldsph to Ldft

            rLdsph = dl_se_ft(rDist, rhte, rhre, rap, rFreq, omega, rpol);

            return rLdsph;

        } else {
            // calculate the smallest clearance between the curved-Earth path and
            // the ray between the antennas, hse

            double rc = (rhte - rhre) / (rhte + rhre);        // Eq (25d)
            double rm = 250 * rDist * rDist / (rap * (rhte + rhre));        // eq (25e)

            double rb = 2 * sqrt((rm + 1) / (3 * rm)) * cos(PI / 3.0 + 1.0 / 3.0 * acos(3 * rc / 2 * sqrt(3 * rm / pow(rm + 1, 3.0))));   // Eq (25c)

            double rdse1 = rDist / 2 * (1 + rb);           // Eq (25a)
            double rdse2 = rDist - rdse1;            // Eq (25b)

            double rhse = (rhte - 500 * rdse1 * rdse1 / rap) * rdse2 + (rhre - 500 * rdse2 * rdse2 / rap) * rdse1;
            rhse = rhse / rDist;                // Eq (24)

            // Calculate the required clearance for zero diffraction loss

            double rhreq = 17.456 * sqrt(rdse1 * rdse2 * rLambda / rDist);     // Eq (26)

            if (rhse > rhreq) {
                rLdsph = 0.0;
                return rLdsph;
            } else {
                // calculate the modified effective Earth radius aem, which gives
                // marginal LoS at distance rDist

                double raem = 500 * pow(rDist / (sqrt(rhte) + sqrt(rhre)), 2.0);     // Eq (27)

                // Use the method in Sec. 4.2.2.1 for adft ) aem to obtain Ldft

                double rLdft = dl_se_ft(rDist, rhte, rhre, raem, rFreq, omega, rpol);

                if (rLdft < 0) {
                    rLdsph = 0.0;
                    return rLdsph;
                }
                else {
                    rLdsph = (1 - rhse / rhreq) * rLdft;     // Eq (28)
                }
            }
        }

        return rLdsph;
    }


double dl_bull(vector<double> & rDisti, vector<double>  & rhi, double rhts, double rhrs, double rAe, double rFreq){
        /**
         This function computes the Bullington part of the diffraction loss
         as defined in ITU-R P.452-17 in 4.2.1

         Input parameters:
         rDisti       -   vector of distances di of the i-th profile point (km)
         rhi          -   vector of heights hi of the i-th profile point (meters above mean sea level. Both vectors contain n+1 profile points
         rhts         -   transmitter antenna height in meters above sea level (i=0)
         rhrs         -   receiver antenna height in meters above sea level (i=n)
         rAe          -   the effective earth radius in kilometers
         rFreq        -   frequency expressed in GHz

         Output parameters:
         Lbull   -   Bullington diffraction loss for a given path

         Rev   Date        Author                          Description
         -------------------------------------------------------------------------------
         v0    31MAR22     Ivica Stevanovic, OFCOM         First implementation in C++
         */

        double rHi = 0;                //
        double rLambda = 0.2998/ rFreq;  //rFreq should be in GHz
        double rCe = 1.0/rAe;          // Effective Earth curvature
        double rLuc = 0;
        int n = rDisti.size();
        double rDist = rDisti[n-1]-rDisti[0];

        // Find the intermediate profile point with the highest slope of the line
        // from the transmitter to the point

        double rStim = (rhi[1] + 500.0*rCe*rDisti[1]*(rDist - rDisti[1]) - rhts)/rDisti[1];  // Eq (14)

        // we assume flat terrain here rHi = 0.0

        for (int i = 2; i < rDisti.size(); i++){

            double rStimi = ((rhi[i] + (500.0 * rCe * rDisti[i] * (rDist - rDisti[i]))) - rhts) / rDisti[i];  // Eq (14)

            rStim = max(rStim, rStimi);
        }

        // Calculate the slope of the line from transmitter to receiver assuming a LoS path

        double rStr = (rhrs - rhts)/rDist;                                         // Eq (15)

        if (rStim < rStr) { // Case 1, Path is LoS

            // Find the intermediate profile point with the highest diffraction  parameter nu:
            double rnumax = rhi[1] + 500 * rCe * rDisti[1] * (rDist - rDisti[1]) - (rhts * (rDist - rDisti[1]) + rhrs * rDisti[1]) / rDist;
            rnumax = rnumax * sqrt(0.002 * rDist / (rLambda * rDisti[1] * (rDist - rDisti[1])));

            for (int i = 2; i < rDisti.size(); i++) {

                double rnumai = rhi[i] + 500 * rCe * rDisti[i] * (rDist - rDisti[i]) - (rhts * (rDist - rDisti[i]) + rhrs * rDisti[i]) / rDist;
                rnumai = rnumai * sqrt(0.002 * rDist / (rLambda * rDisti[i] * (rDist - rDisti[i])));

                rnumax = max(rnumai, rnumax);
            }

            rLuc = 0.0;

            if (rnumax > -0.78) {
                rLuc = knife_edge_loss(rnumax);   // Eq (13), (17)
            }

        }

        else
        {

            // Path is transhorizon

            // Find the intermediate profile point with the highest slope of the  line from the receiver to the point
            double rSrim = (rhi[1] + 500.0 * rCe * rDisti[1] * (rDist - rDisti[1]) - rhrs) / (rDist - rDisti[1]);  // Eq (18)

            for (int i = 2; i < rDisti.size(); i++) {

                double rSrimi = (rhi[i] + 500.0 * rCe * rDisti[i] * (rDist - rDisti[i]) - rhrs) / (rDist - rDisti[i]);  // Eq (18)

                rSrim = max(rSrim, rSrimi);
            }

            // Calculate the distance of the Bullington point from the transmitter:

            double rdbp = (rhrs - rhts + rSrim * rDist) / (rStim + rSrim);                // Eq (19)

            // Calculate the diffraction parameter, nub, for the Bullington point

            double rnub = (rhts + rStim * rdbp - (rhts * (rDist - rdbp) + rhrs * rdbp) / rDist);
            rnub = rnub * sqrt(0.002 * rDist / (rLambda * rdbp * (rDist - rdbp)));    // Eq (20)

            // The knife-edge loss for the Bullington point is given by

            rLuc = 0.0;
            if (rnub > -0.78) {
                rLuc = knife_edge_loss(rnub);   // Eq (13), (21)
            }

        }

        // For Luc calculated using either (17) or (21), Bullington diffraction loss
        // for the path is given by

        return rLuc + (1.0 - exp(-rLuc/6.0))*(10+0.02*rDist);         // Eq (22)



    }
    
double dl_delta_bull( vector<double> & d, vector<double> & h, double hts, double hrs, double hstd, double hsrd, double ap, double f, double omega, int pol ) {
        //dl_delta_bull Complete 'delta-Bullington' diffraction loss model P.452-17
        //   function Ld = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega,pol )
        //
        //   This function computes the complete 'delta-Bullington' diffraction loss
        //   as defined in ITU-R P.452-17 (Section 4.5.3)
        //
        //     Input parameters:
        //     d       -   vector of distances di of the i-th profile point (km)
        //     h       -   vector of heights hi of the i-th profile point (meters
        //     above mean sea level. Both vectors contain n+1 profile points
        //     hts     -   transmitter antenna height in meters above sea level (i=0)
        //     hrs     -   receiver antenna height in meters above sea level (i=n)
        //     hstd    -   Effective height of interfering antenna (m amsl) c.f. 5.1.6.3
        //     hsrd    -   Effective height of interfered-with antenna (m amsl) c.f. 5.1.6.3
        //     ap      -   the effective Earth radius in kilometers
        //     f       -   frequency expressed in GHz
        //     omega   -   the fraction of the path over sea
        //     pol     -   polarization (1 - horizontal, 2 - vertical)
        //
        //     Output parameters:
        //     Ld     -   diffraction loss for the general patha according to
        //                Section 4.2.3 of ITU-R P.452-17.

        //
        //     Example:
        //     Ld = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega, pol )
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         Initial version in C++


        // Use the method in 4.2.1 for the actual terrain profile and antenna
        // heights. Set the resulting Bullington diffraction loss for the actual
        // path to Lbulla

        double Lbulla = dl_bull(d, h, hts, hrs, ap, f);

        // Use the method in 4.2.1 for a second time, with all profile heights hi
        // set to zero and modified antenna heights given by

        double hts1 = hts - hstd;   // eq (38a)
        double hrs1 = hrs - hsrd;   // eq (38b)
        int n = d.size();

        vector<double> h1;
        h1.reserve(n);

        for (int i = 0; i < n; i++) {
            h1.push_back( 0.0 );
        }

        // where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
        // resulting Bullington diffraction loss for this smooth path to Lbulls

        double Lbulls = dl_bull(d, h1, hts1, hrs1, ap, f);

        // Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
        // for the actual path length (dtot) with

        double hte = hts1;             // eq (39a)
        double hre = hrs1;             // eq (39b)
        double dtot = d[n - 1] - d[0];

        double Ldsph = dl_se(dtot, hte, hre, ap, f, omega, pol);
        
        // Diffraction loss for the general paht is now given by

        return Lbulla + max(Ldsph - Lbulls, 0.0);  // eq (40)
    }


void dl_p( vector<double> & d, vector<double> & h, double hts, double hrs, double hstd, double hsrd, double f, double omega, double p, double b0, double DN, int pol, 
              double & Ldp, double & Ld50 ) {
        //dl_p Diffraction loss model not exceeded for p// of time according to P.452-17
        //   function [Ldp, Ld50] = dl_p( d, h, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN )
        //
        //   This function computes the diffraction loss not exceeded for p// of time
        //   as defined in ITU-R P.452-17 (Section 4.5.4)
        //
        //     Input parameters:
        //     d       -   vector of distances di of the i-th profile point (km)
        //     h       -   vector of heights hi of the i-th profile point (meters
        //                 above mean sea level. Both vectors contain n+1 profile points
        //     hts     -   transmitter antenna height in meters above sea level (i=0)
        //     hrs     -   receiver antenna height in meters above sea level (i=n)
        //     hstd    -   Effective height of interfering antenna (m amsl) c.f. 5.1.6.3
        //     hsrd    -   Effective height of interfered-with antenna (m amsl) c.f. 5.1.6.3
        //     f       -   frequency expressed in GHz
        //     omega   -   the fraction of the path over sea
        //     p       -   percentage of time
        //     b0      -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
        //     DN      -   the average radio-refractive index lapse-rate through the
        //                 lowest 1 km of the atmosphere. Note that DN is positive
        //                 quantity in this procedure
        //     pol     -   polarization (1-horizontal, 2-vertical)
        //
        //     Output parameters:
        //     Ldp    -   diffraction loss for the general path not exceeded for p // of the time
        //                according to Section 4.2.4 of ITU-R P.452-17.
        //                Ldp(1) is for the horizontal polarization
        //                Ldp(2) is for the vertical polarization
        //     Ld50   -   diffraction loss for p = 50//
        //
        //     Example:
        //     [Ldp, Ld50] = dl_p( d, h, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN )
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         Initial version in C++


        // Use the method in 4.2.3 to calculate diffraction loss Ld for effective
        // Earth radius ap = ae as given by equation (6a). Set median diffraction
        // loss to Ldp50
        double ae, ab;
        earth_rad_eff(DN, ae, ab);

        double ap = ae;

        Ld50 = dl_delta_bull(d, h, hts, hrs, hstd, hsrd, ap, f, omega, pol);

        Ldp = 0.0;

        if (p == 50) {
            Ldp = Ld50;

        }


        if (p < 50) {

            // Use the method in 4.2.3 to calculate diffraction loss Ld for effective
            // Earth radius ap = abeta, as given in equation (6b). Set diffraction loss
            // not exceeded for beta0// time Ldb = Ld

            ap = ab;

            double Ldb;
            Ldb = dl_delta_bull(d, h, hts, hrs, hstd, hsrd, ap, f, omega, pol);

            // Compute the interpolation factor Fi
            double Fi = 1.0;
            if (p > b0) {

                Fi = inv_cum_norm(p / 100) / inv_cum_norm(b0 / 100.0);   // eq (41a)
            }


            // The diffraction loss Ldp not exceeded for p// of time is now given by

            Ldp = Ld50 + Fi * (Ldb - Ld50);   // eq (42)

        }

    }


void closs_corr(double f, vector<double> & d, vector <double> & h, vector<int> & zone, double htg, double hrg, double ha_t, double ha_r, double dk_t, double dk_r,
                int & index1, int & index2, double & htgc, double & hrgc, double & Aht, double & Ahr) {
        // [index1, index2, htgc, hrgc, Aht, Ahr] = closs_corr(f, d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r)
        //closs clutter loss correction according to P.452-17
        //
        //   This function computes the height-gain correction as defined in ITU-R P.452-17 (Section 4.5.4)
        //
        //     Input parameters:
        //     f       -   Frequency (GHz)
        //     d       -   vector of distances di of the i-th profile point (km)
        //     h       -   vector of heights hi of the i-th profile point (meters
        //                 above mean sea level. Both vectors contain n+1 profile points
        //     zone    -   Zone type: Coastal land (1), Inland (2) or Sea (3)
        //     htg     -   Tx Antenna center heigth above ground level (m)
        //     hrg     -   Rx Antenna center heigth above ground level (m)
        //     ha_t    -   Nominal clutter height at the transmitting end (m, agl)
        //     ha_r    -   Nominal clutter height at the receiving end (m, agl)
        //     dk_t    -   distance from nominal clutter point to the Tx antenna (km)
        //     dk_r    -   distance from nominal clutter point to the Rx antenna (km)
        //
        //     Output parameters:
        //     index1  -
        //     index2  -   Indices defining the profile points between terminal clutter   
        //     htgc    -   Tx Antenna center heigth above ground level (m) in the height-gain model
        //     hrgc    -   Rx Antenna center heigth above ground level (m) in the height-gain model
        //     Aht     -   Additional losses to account for clutter shielding the
        //     Ahr         transmitter and receiver. These should be set to zero if there is no
        //                 such shielding
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         Initial version in C++


        index1 = 0;
        index2 = d.size()-1;

        htgc = htg;
        hrgc = hrg;

        Aht = 0;
        Ahr = 0;

        double ha = ha_t;
        double dk = dk_t;

        if (ha > htg) {



            double Ffc = 0.25 + 0.375 * (1 + tanh(7.5 * (f - 0.5)));  // (57a)

            Aht = 10.25 * Ffc * exp(-dk) * (1 - tanh(6 * (htg / ha - 0.625))) - 0.33; // (57)

            for (int kk=0; kk< d.size(); kk++) {

                if (d[kk] >= dk) {
                    index1 = kk;
                    break;
                }

            }

            htgc = ha_t;

        }

        ha = ha_r;

        dk = dk_r;

        if (ha > hrg) {

            double Ffc = 0.25 + 0.375 * (1 + tanh(7.5 * (f - 0.5)));  // (57a)

            Ahr = 10.25 * Ffc * exp(-dk) * (1 - tanh(6 * (hrg / ha - 0.625))) - 0.33;  // (57)

            for (int kk = d.size()-1; kk >=0; kk--){
                if (d[kk] <= d[d.size()-1] - dk ){
                    index2 = kk;
                    break;
                }
            }

            hrgc = ha_r;
        }

        // make sure there are at least three terrain profile points between Tx and Rx clutter sites
        if (index2 - index1  < 4) {
            
            throw std::range_error("There must be at least three profile points between Tx and Rx terminal clutter locations.");
            
        }
    }


    bool check_limit(double var, double low, double hi, std::string name) {
        if ((var < low) || (var > hi)){ 
            std::stringstream message;
            message << name << " = " << var <<  " is outside the limits: [" << low <<  ", " << hi << "]" << std::endl;
            throw std::invalid_argument(message.str());
            return false;
        }
        return true;
    }


    bool check_value(vector<int> & var, vector<int> & dvals, std::string name) {
        int n = var.size();
        int nd = dvals.size();
        bool res = false;
        
        for (int i = 0; i < n; i++){
            res = false;
            for (int id = 0; id < nd; id++){
                res = res || (var[i] == dvals[id]) ;
            }
            if (res == false){
                std::stringstream message;
                std::ostringstream oss;
                std::copy(dvals.begin(), dvals.end(), std::ostream_iterator<int>(oss, " "));
                
                message << name << " may only contain the following values: " << oss.str() << std::endl;
                
                throw std::invalid_argument(message.str());
                return false;
            }

        }
        return true;
    }

    bool check_value(int var, vector<int> & dvals, std::string name) {
        int nd = dvals.size();
        bool res = false;
        
        res = false;
        for (int id = 0; id < nd; id++){
            res = res || (var == dvals[id]) ;
        }
        if (res == false){
            std::stringstream message;
            std::ostringstream oss;
            std::copy(dvals.begin(), dvals.end(), std::ostream_iterator<int>(oss, " "));
            
            message << name << " may only contain the following values: " << oss.str() << std::endl;
            
            throw std::invalid_argument(message.str());
            return false;
        }


        return true;
    }

    void verify_inputs(double f, double p, vector<double> & d, vector<double> & h, vector<double> &g, vector<int> & zone, double htg, double hrg, double phit_e, double phit_n, double phir_e, double phir_n,  double Gt, double Gr, int pol, double dct, double dcr, double press, double temp) {

        // values allowed for polarization
        vector<int> pol_ref;
        pol_ref.push_back(1);
        pol_ref.push_back(2);

        
        // values allowed for radio-climatic zone
        vector<int> zone_ref;
        zone_ref.push_back(1);
        zone_ref.push_back(2);
        zone_ref.push_back(3);

        // check the limits

        check_limit(f, 0.1, 50.0, "f [GHz]");
        check_limit(p, 0.001, 50, "p [%]");
        check_limit(phit_n, -89.9999, 89.9999, "phit_n [deg]");
        check_limit(phir_n, -89.9999, 89.9999, "phir_n [deg]");
        check_limit(phit_e, -180, 360, "phit_e [deg]");
        check_limit(phir_e, -180, 360, "phir_e [deg]");
        check_limit(dct, 0, std::numeric_limits<double>::max() , "dct [km]");
        check_limit(dcr, 0, std::numeric_limits<double>::max() , "dcr [km]");
        check_value(pol, pol_ref, "Polarization (pol) ");
        check_value(zone, zone_ref, "Radio-climatic zone (zone) ");
        
        if (d[0] > 0.0) {
            
            throw std::invalid_argument("d[0] must be equal to zero. ");
            
        }
    }

void tropospheric_path(double dt, double hts, double hrs, double theta_e, double theta_tpos, double theta_rpos, double phi_re, double phi_te, double phi_rn, double phi_tn, double Re, double & d_tcv, double & phi_cve, double & phi_cvn) {
        //trophospheric path segments according to ITU-R P.2001-5
        // This function computes tropospheric path segments as described in Section
        // 3.9 of Recommendation ITU-R P.2001-5
        //
        // Input parameters:
        // dt        -   Path length (km)
        // hts, hrs  -   Tx/Rx antenna heights above means sea level (m)
        // theta_e   -   Angle subtended by d km at the center of a sphere of effective earth radius (rad)
        // theta_tpos-   Interfering antenna horizon elevation angle limited to be positive (mrad)
        // theta_rpos-   Interfered-with antenna horizon elevation angle limited to be positive (mrad)
        //               hts = htg + h[0]
        // phi_re    -   Receiver longitude, positive to east (deg)
        // phi_te    -   Transmitter longitude, positive to east (deg)
        // phi_rn    -   Receiver latitude, positive to north (deg)
        // phi_tn    -   Transmitter latitude, positive to north (deg)
        // Re        -   Average Earth radius (km)
        //
        // Output parameters:
        // d_tcv     -   Horizontal path length from transmitter to common volume (km)
        // phi_cve   -   Longitude of the common volume
        // phi_cvn   -   Latitude of the common volume
        //
        // Rev   Date        Author                          Description
        // -------------------------------------------------------------------------------
        // v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version
       
        // Horizontal path lenght from transmitter to common volumne (3.9.1a)

        d_tcv = (dt * tan(0.001 * theta_rpos + 0.5 * theta_e) - 0.001 * (hts - hrs)) /
                (tan(0.001 * theta_tpos + 0.5 * theta_e) + tan(0.001 * theta_rpos + 0.5 * theta_e));

        // Limit d_tcv such that 0 <= dtcv <= dt

        if (d_tcv < 0) {
            d_tcv = 0;
        }
        if (d_tcv > dt) {
            d_tcv = dt;
        }


        // Calculate the longitude and latitude of the common volumne from the
        // transmitter and receiver longitudes and latitudes using the great circle
        // path method of Attachment H by seting d_pnt = d_tcv
        double b2r, dgc;

        great_circle_path(phi_re, phi_te, phi_rn, phi_tn, Re, d_tcv, phi_cve, phi_cvn, b2r, dgc);

    }

void  surface_altitude_cv(vector<double> h, vector <double> d, double d_tcv, double & hs) {
        // surface_altitude_cv altitude on the surface of the Earth below common volume
        //     surface_altitude_cv(h, d, d_tcv, hs)
        //     This function computes the altitude of the point at the surface of
        //     the Earth below common volume
        //
        //     Input arguments:
        //           d       -   vector of distances in the path profile (km)
        //           h       -   vector of heights (masl)
        //           d_ctv   -   horizontal path length from transmitter to common volume computed using (3.9.1a)
        //
        //     Output arguments:
        //           hs      -   altitude on the surface of the Earth below common volume (masl)
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    20JUL23     Ivica Stevanovic, OFCOM         First implementation 


        int n = d.size();

        int ii = 0;
        double dmin = abs(d[0] - d_tcv);

        for (int i = 1; i < n; i++) {
            double dmin_trial = abs(d[i] - d_tcv);
            if (dmin_trial < dmin) {
                dmin = dmin_trial;
                ii = i;
            }
        }

        int i1 = 0;
        int i2 = 0;
        hs = 0;

        if (d[ii] == d_tcv) {
            i1 = ii;
            hs = h[i1];
            return ;
        }

        if (d[ii] < d_tcv) {
            i1 = ii;
            i2 = ii + 1;
        } else {
            i2 = ii;
            i1 = ii - 1;
        }

        // apply linear interpolation

        hs = h[i1] + (h[i2] - h[i1]) * (d_tcv - d[i1]) / (d[i2] - d[i1]);
        return;
    }

void tl_troposcatter(P452DigitalMaps maps, double f, double dt, double hts, double hrs, double ae, double thetae, double thetat, double thetar,  double phicvn, double phicve, double Gt, double Gr, double p, double hs, double & Lbs, double & theta) {
        //tl_troposcatter Troposcatter basic transmission loss
        //   This function computes the troposcatter basic transmission loss
        //   as defined in Section 4.3
        //
        //     Input parameters:
        //     maps    -   Object containing Digital Maps
        //     f       -   Frequency GHz
        //     dt      -   Total distance (km)
        //     hts,hrs -   Altitudes of transmitting antenna and receiving antennas in m
        //     ae      -   Effective Earth radius (km)
        //     thetae  -   Angle subtended by d km at centre of spherical Earth (rad)
        //     thetat  -   Tx horizon elevation angle relative to the local horizontal (mrad)
        //     thetar  -   Rx horizon elevation angle relative to the local horizontal (mrad)
        //     phicvn  -   Troposcatter common volume latitude (deg)
        //     phicve  -   Troposcatter common volume longitude (deg)
        //     Gt, Gr  -   Gain of transmitting and receiving antenna in the azimuthal direction
        //                 of the path towards the other antenna and at the elevation angle
        //                 above the local horizontal of the other antenna in the case of a LoS
        //                 path, otherwise of the antenna's radio horizon, for median effective
        //                 Earth radius.

        //     p       -   Percentage of average year for which predicted basic loss
        //                 is not exceeded (%)
        //     hs      -   Height of the Earth's surface above sea level (km)
        //
        //     Output parameters:
        //     Lbs    -   Troposcatter basic transmission loss (dB)
        //     theta  -   Scatter angle (mrad)
        //
        //
        //     Example:
        //     [Lbs, theta] = tl_troposcatter_pdr(f, dt, hts, hrs, ae, the, thetat, thetar, phicvn, phicve,  Gt, Gr, p, hs, Lbs, theta)

        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version


        // Attachment E: Troposcatter

        double fMHz = f*1000;

        // E.2 Climatic classification

        // Find average annual sea-level surface refractivity N0 and radio-refractivity lapse rate dN
        // for the common volume of the link in question using the digital maps at phicve (lon),
        // phicvn (lat) - as a bilinear interpolation
        double dN = maps.GetDN50(phicve, phicvn);
        double N0 = maps.GetN050(phicve, phicvn);


        // E.3 Calculation of tropocscatter basic transmission loss
        // Step2: Calculate the scatter angle theta

        theta = 1000*thetae + thetat + thetar;     // mrad    (145)


        // Step 3: Estimate the aperture-to-median coupling loss Lc (11)

        double Lc = 0.07 * exp(0.055* (Gt + Gr));    // dB    (45a)


        // Step 4: Estimate the average annual transmission loss associated with
        // troposcatter not exceeded for p% of time (45):

        double hb = 7.35;  //km  scale height set to the global mean

        double beta = dt/(2.0*ae) + thetar/1000.0 + (hrs-hts)/(1000.0*dt);  //(45e)

        double h0 = hts/1000.0 + dt*sin(beta)/(sin(theta/1000.0)) *(0.5* dt * sin(beta)/(ae*sin(theta/1000.0))+sin(thetat/1000.0));   //(45d)

        double Yp = 0.035*N0*exp(-h0/hb)*pow( (-log10(p/50.0)), (0.67));

        if (p >= 50) {

            Yp = -0.035 * N0 * exp(-h0 / hb) * pow( (-log10((100.0 - p) / 50.0)),  (0.67));

        }

        double F = 0.18*N0*exp(-hs/hb) - 0.23*dN;   //(45b)

        Lbs = F + 22.0*log10(fMHz) + 35.0*log10(theta) + 17.0*log10(dt) + Lc - Yp;    // (45)

        return;

    }



    double tl_p452(P452DigitalMaps maps, double f, double p, vector<double> & d, vector<double> & h, vector<double> & g, vector<int> & zone, double htg, double hrg, double phit_e, double phit_n, double phir_e, double phir_n,  double Gt, double Gr, int pol, double dct, double dcr, double press, double temp) {
        //tl_p452 basic transmission loss according to P.452-18
        //  Lb = tl_p452(maps, f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, diffraction_only );
        //
        //   This is the MAIN function that computes the basic transmission loss not exceeded for p% of time
        //   as defined in ITU-R P.452-18 (clear air portion).
        //
        //     Input parameters:
        //     maps    -   Object containing all the Digital Maps (DN50, N050) necessary for computation
        //     f       -   Frequency (GHz)
        //     p       -   Required time percentage for which the calculated basic
        //                 transmission loss is not exceeded
        //     d       -   vector of distances di of the i-th profile point (km)
        //     h       -   vector of heights hi of the i-th profile point (meters
        //                 above mean sea level. Both vectors contain n+1 profile points
        //     g       -   vector of clutter + terrain profile heights gi along the path gi = hi + Ri (masl)
        //                 where Ri is the (representative clutter height)
        //     zone    -   Zone type: Coastal land (1), Inland (2) or Sea (3)
        //     htg     -   Tx Antenna center heigth above ground level (m)
        //     hrg     -   Rx Antenna center heigth above ground level (m)
        //     phit_e  -   Tx Longitude (degrees)
        //     phit_n  -   Tx Latitude  (degrees)
        //     phir_e  -   Rx Longitude (degrees)
        //     phir_n  -   Rx Latitude  (degrees)
        //     Gt, Gr  -   Antenna gain in the direction of the horizon along the
        //                 great-circle interference path (dBi)
        //     pol     -   polarization of the signal (1) horizontal, (2) vertical
        //     dct     -   Distance over land from the transmit and receive
        //     dcr         antennas to the coast along the great-circle interference path (km).
        //                 Set to zero for a terminal on a ship or sea platform
        //     press   -   Dry air pressure (hPa)
        //     temp    -   Air temperature (degrees C)
        //     diffraction_only - computes the basic transmission loss from diffraction only (boolean)
        //
        //     Output parameters:
        //     Lb     -   basic  transmission loss according to P.452-18
        //
        //     Example:
        //     Lb = tl_p452(f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, diffraction_only );

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    31MAR22     Ivica Stevanovic, OFCOM         Initial C++ implementation
        //     v1    20JUL23     Ivica Stevanovic, OFCOM         Aligned with Rec. ITU-R P.452-18 (distributed clutter model, harmonized troposcatter model)


        // Redistribution and use in source and binary forms, with or without
        // modification, are permitted provided that the following conditions are
        // met:
        //
        //     * Redistributions of source code must retain the above copyright
        //       notice, this list of conditions and the following disclaimer.
        //     * Redistributions in binary form must reproduce the above copyright
        //       notice, this list of conditions and the following disclaimer in
        //       the documentation and/or other materials provided with the distribution
        //
        //
        ////
        // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
        // ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
        // LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
        // CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
        // SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
        // INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
        // CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
        // ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
        // POSSIBILITY OF SUCH DAMAGE.
        //
        // THE AUTHORS AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
        ////

        // verify input arguments values and limits
        try{
            verify_inputs(f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp);
        }
        
        catch (const std::exception& e) {
            
            std::cout << e.what(); 
            return -1000.0;
        }   

        int n = d.size();
        int ks = 0;
        int ke = n;

        // Apply the condition in Step 4: Radio profile 
        // gi is the terrain height in metres above sea level for all the points at a distance from transmitter or receiver less than 50 m.

        for (int k = 0; k < n; k++) {
            if (d[k] < 50.0 / 1000.0) {
                g[k] = h[k];
            } else {
                break;
            }
        }

        double endVal = d[n-1] - 50.0/1000.0;
        for (int k = n - 1; k >=0 ; k--) {
            if ( d[k] > endVal) {
                g[k] = h[k];
            } else {
                break;
            }
        }


        // Compute the path profile parameters
    
        // Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
        
        int zone_r = 12;
        double dtm = longest_cont_dist(d, zone, zone_r);

        // Compute  dlm     -   the longest continuous inland section of the great-circle path (km)
        zone_r = 2;
        double dlm = longest_cont_dist(d, zone, zone_r);

      // Calculate the longitude and latitude of the mid-point of the path, phim_e,
        // and phim_n for dpnt = 0.5dt
        double Re = 6371;
        double dtot = d[n - 1] - d[0];
        double dpnt = 0.5*dtot;
        double phim_e, phim_n, bt2r, dgc;

        great_circle_path(phir_e, phit_e, phir_n, phit_n, Re, dpnt, phim_e, phim_n, bt2r, dgc);

        // Find radio-refractivity lapse rate dN
        // using the digital maps at phim_e (lon), phim_n (lat) - as a bilinear interpolation

        double DN = maps.GetDN50(phim_e, phim_n);
        double N0 = maps.GetN050(phim_e, phim_n);

        // Compute b0
        double b0 = beta0(phim_n, dtm, dlm);

        double ae, ab;
        earth_rad_eff(DN, ae, ab);

        
        // Compute the path fraction over see

        double omega = path_fraction_sea(d, zone, 3);

        double hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta;
        int pathtype;
        
        smooth_earth_heights(d, h, htg, hrg, ae, f, hst, hsr, hstd, hsrd, hte, hre, hm, dlt,dlr, theta_t, theta_r, theta, pathtype);
        
        //Tx and Rx antenna heights above mean sea level amsl (m)
        
        double hts = h[0] + htg;
        double hrs = h[n - 1] + hrg;

        // Effective Earth curvature Ce (km^-1)

        double Ce = 1.0 / ae;

        // Wavelength in meters

        double lambda = 0.2998 / f;


        // Find the intermediate profile point with the highest slope of the line
        // from the transmitter to the point

        double Stim = ((h[1] + 500 * Ce * d[1] * (dtot - d[1]) - hts) / d[1]);

        for (int i = 2; i < n - 1; i++) {
            Stim = std::max(Stim, (h[i] + 500 * Ce * d[i] * (dtot - d[i]) - hts) / d[i]);       // Use hi instead of gi in Eq (14)
        }

        // Calculate the slope of the line from transmitter to receiver assuming a
        // LoS path

        double Str = (hrs - hts) / dtot;                                         // Eq (15)

        // Calculate an interpolation factor Fj to take account of the path angular
        // distance (58)

        double THETA = 0.3;
        double KSI = 0.8;

        // changed the definition for Fj on 15DEC16.
        //Fj = 1.0 - 0.5*( 1.0 + tanh(3.0 * KSI * (theta-THETA)/THETA) );
        double Fj = 1.0 - 0.5 * (1.0 + tanh(3.0 * KSI * (Stim - Str) / THETA));

        // Calculate an interpolation factor, Fk, to take account of the great
        // circle path distance:

        double dsw = 20;
        double kappa = 0.5;

        double Fk = 1.0 - 0.5 * (1.0 + tanh(3.0 * kappa * (dtot - dsw) / dsw)); // eq (59)

        double Lbfsg, Lb0p, Lb0b;

        double d3D = sqrt(dtot*dtot + pow( (hts-hrs)/1000.0 , 2.0 ) );
        pl_los(d3D, f, p, b0, omega, temp, press, dlt, dlr, Lbfsg, Lb0p, Lb0b);

        double Ldp, Ld50;
        dl_p(d, g, hts, hrs, hstd, hsrd, f, omega, p, b0, DN, pol, Ldp, Ld50);

        // The median basic transmission loss associated with diffraction Eq (43)

        double Lbd50 = Lbfsg + Ld50;

        // The basic tranmission loss associated with diffraction not exceeded for
        // p% time Eq (44)

        double Lbd = Lb0p + Ldp;

        // A notional minimum basic transmission loss associated with LoS
        // propagation and over-sea sub-path diffraction

        double Lminb0p = Lb0p + (1 - omega) * Ldp;

        double Fi = 1;

        if (p >= b0) {

            Fi = inv_cum_norm(p / 100.0) / inv_cum_norm(b0 / 100.0);   //eq (41a)

            Lminb0p = Lbd50 + (Lb0b + (1 - omega) * Ldp - Lbd50) * Fi;   // eq (60)

        }

        // Calculate a notional minimum basic transmission loss associated with LoS
        // and transhorizon signal enhancements

        double eta = 2.5;

        double Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, ae, b0);

        double Lminbap = eta * log(exp(Lba / eta) + exp(Lb0p / eta));    // eq (61)

        // Calculate a notional basic transmission loss associated with diffraction
        // and LoS or ducting/layer reflection enhancements

        double Lbda = Lbd;

        if (Lminbap <= Lbd) {
            Lbda = Lminbap + (Lbd - Lminbap) * Fk;
        }

        // Calculate a modified basic transmission loss, which takes diffraction and
        // LoS or ducting/layer-reflection enhancements into account

        double Lbam = Lbda + (Lminb0p - Lbda) * Fj;   // eq (63)

       // Calculate the basic transmission loss due to troposcatter not exceeded
        // for any time percentage p

        /* The path length expressed as the angle subtended by d km at the center of
        // a sphere of effective Earth radius ITU-R P.2001-4 (3.5.4)

        double theta_e = dtot/ae; // radians

        //Calculate the horizon elevation angles limited such that they are positive

        double theta_tpos = std::max(theta_t, 0.0);                   // Eq (3.7.11a) ITU-R P.2001-5
        double theta_rpos = std::max(theta_r, 0.0);                   // Eq (3.7.11b) ITU-R P.2001-5

        double d_tcv, phi_cve, phi_cvn;

        tropospheric_path(dtot, hts, hrs, theta_e, theta_tpos, theta_rpos, phir_e, phit_e, phir_n, phit_n, Re, d_tcv, phi_cve, phi_cvn);
       //     std::cout << "d_tcv = " << d_tcv << "    phi_cve = " << phi_cve  <<"      phi_cvn = " << phi_cvn << std::endl;
        double Hs;
        surface_altitude_cv(h, d, d_tcv, Hs);
        Hs = Hs/1000.0;

        double Lbs, Thetas;
        tl_troposcatter(maps, f, dtot, hts, hrs, ae, theta_e, theta_t, theta_r, phi_cvn, phi_cve, Gt, Gr, p, Hs, Lbs, Thetas);

        // To avoid under-estimating troposcatter for short paths, limit Lbs as in ITU-R P.2001-5 (E.17)

        Lbs = std::max(Lbs, Lbfsg);
        */

        double Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr);

        // Calculate the final transmission loss not exceeded for p% time


        return -5.0 * log10(pow(10.0, -0.2 * Lbs) + pow(10.0, -0.2 * Lbam));  // eq (64)

     
    }    
