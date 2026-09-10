// validate_p452.cpp
//
// C++ port of validate_p452.m (see
// https://github.com/eeveetza/p452/blob/main/matlab/validate_p452.m).
//
// Unlike the previous test1.cpp .. test17.cpp files (produced once, ahead
// of time, by tests/create_cpp_validation_p452.m and then checked in with
// hardcoded numbers), this test reads the very same .csv files the MATLAB
// script uses - tests/validation_examples/profiles/*.csv and the matching
// tests/validation_examples/results/*.csv - at test run time. It computes
// both the intermediate path-profile parameters (ae, dtot, hts, hrs, ...)
// and the final/intermediate transmission-loss terms (Lbfsg, Lb0p, Lb0b,
// Ldsph, Ld50, Ldp, Lbs, Lba, Lb) exactly like validate_p452.m does, and
// compares them against the reference values in the .csv files.
//
// Whenever the MATLAB repository's validation_examples/ folder is
// updated, this test picks up the change automatically the next time it
// runs - no regeneration step is needed.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "../include/P452.hpp"
#include "csv_utils.hpp"

namespace fs = std::filesystem;

namespace {

constexpr double kTol = 1e-6;

// Walks up from the current working directory until it finds relPath,
// the same strategy the previous hand-generated tests used to locate
// src/maps/ regardless of where ctest happens to set the working
// directory.
std::string resolveRepoPath(const std::string& relPath) {
  auto baseDir = fs::current_path();
  while (baseDir.has_relative_path()) {
    auto candidate = baseDir / relPath;
    if (fs::exists(candidate)) {
      return candidate.string();
    }
    baseDir = baseDir.parent_path();
  }
  throw std::runtime_error("Could not locate '" + relPath +
                           "' in any parent of " + fs::current_path().string());
}

struct ProfileData {
  std::vector<double> d;  // distance (km)
  std::vector<double> h;  // terrain height (m amsl)
  std::vector<double> r;  // representative clutter height (m)
  std::vector<int> zone;  // 1 - Coastal Land, 2 - Inland, 3 - Sea
};

// Mirrors the "read the path profile" section of validate_p452.m:
// columns are d (km), h (m), clutter height (m), zone letter, zone number.
ProfileData readProfile(const std::string& path) {
  ProfileData profile;
  const auto rows = readCsv(path);

  profile.d.reserve(rows.size());
  profile.h.reserve(rows.size());
  profile.r.reserve(rows.size());
  profile.zone.reserve(rows.size());

  for (const auto& row : rows) {
    profile.d.push_back(std::stod(row[0]));
    profile.h.push_back(std::stod(row[1]));
    profile.r.push_back(std::stod(row[2]));
    profile.zone.push_back(std::stoi(row[4]));
  }

  return profile;
}

struct ResultsData {
  // Per-row (f, p) test cases sharing the same profile and fixed inputs.
  std::vector<double> f;
  std::vector<double> p;

  // Fixed scalar inputs (identical for every row of a given file).
  double htg = 0, hrg = 0;
  double phit_e = 0, phit_n = 0, phir_e = 0, phir_n = 0;
  double Gt = 0, Gr = 0;
  int pol = 0;
  double dct = 0, dcr = 0;
  double press = 0, temp = 0;

  // Reference intermediate path-profile parameters (identical for every
  // row - computed once from the profile and the first frequency, exactly
  // as in validate_p452.m / tl_p452).
  double ae_ref = 0, dtot_ref = 0, hts_ref = 0, hrs_ref = 0;
  double theta_t_ref = 0, theta_r_ref = 0, theta_ref = 0;
  double hm_ref = 0, hte_ref = 0, hre_ref = 0, hstd_ref = 0, hsrd_ref = 0;
  double dlt_ref = 0, dlr_ref = 0;
  int pathtype_ref = 0;
  double dtm_ref = 0, dlm_ref = 0, b0_ref = 0, omega_ref = 0;
  double DN_ref = 0, N0_ref = 0;

  // Reference per-row transmission-loss terms.
  std::vector<double> Lb_ref, Lbfsg_ref, Lb0p_ref, Lb0b_ref;
  std::vector<double> Ldsph_ref, Ld50_ref, Ldp_ref, Lbs_ref, Lba_ref;
};

// Column layout (0-based) of validation_examples/results/*.csv, per the
// header row written out by validate_p452.m:
// profile,f,p,htg,hrg,phit_e,phit_n,phir_e,phir_n,Gt,Gr,pol,dct,dcr,press,
// temp,ae,dtot,hts,hrs,theta_t,theta_r,theta,hm,hte,hre,hstd,hsrd,dlt,dlr,
// path,dtm,dlm,b0,omega,DN,N0,Lb,Lbfsg,Lb0p,Lb0b,Ldsph,Ld50,Ldp,Lbs,Lba
ResultsData readResults(const std::string& path) {
  ResultsData r;
  const auto rows = readCsv(path);
  if (rows.empty()) {
    throw std::runtime_error("No data rows in " + path);
  }

  const size_t n = rows.size();
  r.f.resize(n);
  r.p.resize(n);
  r.Lb_ref.resize(n);
  r.Lbfsg_ref.resize(n);
  r.Lb0p_ref.resize(n);
  r.Lb0b_ref.resize(n);
  r.Ldsph_ref.resize(n);
  r.Ld50_ref.resize(n);
  r.Ldp_ref.resize(n);
  r.Lbs_ref.resize(n);
  r.Lba_ref.resize(n);

  for (size_t i = 0; i < n; ++i) {
    const auto& row = rows[i];
    r.f[i] = std::stod(row[1]);
    r.p[i] = std::stod(row[2]);
    r.Lb_ref[i] = std::stod(row[37]);
    r.Lbfsg_ref[i] = std::stod(row[38]);
    r.Lb0p_ref[i] = std::stod(row[39]);
    r.Lb0b_ref[i] = std::stod(row[40]);
    r.Ldsph_ref[i] = std::stod(row[41]);
    r.Ld50_ref[i] = std::stod(row[42]);
    r.Ldp_ref[i] = std::stod(row[43]);
    r.Lbs_ref[i] = std::stod(row[44]);
    r.Lba_ref[i] = std::stod(row[45]);
  }

  const auto& row0 = rows[0];
  r.htg = std::stod(row0[3]);
  r.hrg = std::stod(row0[4]);
  r.phit_e = std::stod(row0[5]);
  r.phit_n = std::stod(row0[6]);
  r.phir_e = std::stod(row0[7]);
  r.phir_n = std::stod(row0[8]);
  r.Gt = std::stod(row0[9]);
  r.Gr = std::stod(row0[10]);
  r.pol = std::stoi(row0[11]);
  r.dct = std::stod(row0[12]);
  r.dcr = std::stod(row0[13]);
  r.press = std::stod(row0[14]);
  r.temp = std::stod(row0[15]);

  r.ae_ref = std::stod(row0[16]);
  r.dtot_ref = std::stod(row0[17]);
  r.hts_ref = std::stod(row0[18]);
  r.hrs_ref = std::stod(row0[19]);
  r.theta_t_ref = std::stod(row0[20]);
  r.theta_r_ref = std::stod(row0[21]);
  r.theta_ref = std::stod(row0[22]);
  r.hm_ref = std::stod(row0[23]);
  r.hte_ref = std::stod(row0[24]);
  r.hre_ref = std::stod(row0[25]);
  r.hstd_ref = std::stod(row0[26]);
  r.hsrd_ref = std::stod(row0[27]);
  r.dlt_ref = std::stod(row0[28]);
  r.dlr_ref = std::stod(row0[29]);

  const std::string pathStr = row0[30];
  r.pathtype_ref = (pathStr.find("Line of Sight") != std::string::npos) ? 1 : 2;

  r.dtm_ref = std::stod(row0[31]);
  r.dlm_ref = std::stod(row0[32]);
  r.b0_ref = std::stod(row0[33]);
  r.omega_ref = std::stod(row0[34]);
  r.DN_ref = std::stod(row0[35]);
  r.N0_ref = std::stod(row0[36]);

  return r;
}

// validate_p452.m derives "test_result<suffix>.csv" from
// "test_profile<suffix>.csv" by stripping the "test_profile" prefix (the
// leading underscore is kept) and prepending "test_result" instead.
std::string resultsFilenameFor(const std::string& profileFilename) {
  const std::string prefix = "test_profile";
  if (profileFilename.rfind(prefix, 0) != 0) {
    throw std::runtime_error("Unexpected profile filename: " + profileFilename);
  }
  return "test_result" + profileFilename.substr(prefix.size());
}

}  // namespace

TEST(P452Validation, AllValidationExamples) {
  const std::string mapsPath = resolveRepoPath("src/maps/");
  P452DigitalMaps maps(mapsPath);

  const std::string profilesDir = resolveRepoPath("tests/validation_examples/profiles/");
  const std::string resultsDir = resolveRepoPath("tests/validation_examples/results/");

  std::vector<std::string> profileFiles;
  for (const auto& entry : fs::directory_iterator(profilesDir)) {
    if (entry.path().extension() == ".csv") {
      profileFiles.push_back(entry.path().filename().string());
    }
  }
  std::sort(profileFiles.begin(), profileFiles.end());

  ASSERT_FALSE(profileFiles.empty())
      << "No validation profiles found in " << profilesDir;

  for (const auto& profileFile : profileFiles) {
    SCOPED_TRACE("Profile file: " + profileFile);

    const ProfileData profile = readProfile(profilesDir + profileFile);
    const ResultsData ref = readResults(resultsDir + resultsFilenameFor(profileFile));

    const size_t n = profile.d.size();
    ASSERT_EQ(n, profile.h.size());
    ASSERT_EQ(n, profile.r.size());
    ASSERT_EQ(n, profile.zone.size());
    ASSERT_FALSE(ref.f.empty());

    // g = terrain + clutter height, with the Step 4 override: within 50 m
    // of either terminal, g falls back to the bare terrain height h.
    std::vector<double> g(n);
    for (size_t i = 0; i < n; ++i) {
      g[i] = profile.h[i] + profile.r[i];
    }
    for (size_t i = 0; i < n && profile.d[i] < 50.0 / 1000.0; ++i) {
      g[i] = profile.h[i];
    }
    const double endVal = profile.d[n - 1] - 50.0 / 1000.0;
    for (size_t i = n; i-- > 0 && profile.d[i] > endVal;) {
      g[i] = profile.h[i];
    }

    // --- Path profile parameters (mirrors validate_p452.m) ---

    std::vector<double> d = profile.d;   // non-const copies: the C++
    std::vector<double> h = profile.h;   // functions below take
    std::vector<int> zone = profile.zone; // non-const references.

    int zoneR = 12;
    const double dtm = longest_cont_dist(d, zone, zoneR);
    zoneR = 2;
    const double dlm = longest_cont_dist(d, zone, zoneR);

    const double Re = 6371.0;
    const double dtot = d[n - 1] - d[0];
    const double dpnt = 0.5 * dtot;
    double phim_e, phim_n, bt2r, dgc;
    great_circle_path(ref.phir_e, ref.phit_e, ref.phir_n, ref.phit_n, Re, dpnt,
                      phim_e, phim_n, bt2r, dgc);

    const double DN = maps.GetDN50(phim_e, phim_n);
    const double N0 = maps.GetN050(phim_e, phim_n);

    const double b0 = beta0(phim_n, dtm, dlm);

    double ae, ab;
    earth_rad_eff(DN, ae, ab);

    // Smooth-Earth heights are computed once, from the first frequency in
    // the file, and reused for every (f, p) row - exactly as
    // validate_p452.m does.
    double hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta;
    int pathtype;
    smooth_earth_heights(d, h, ref.htg, ref.hrg, ae, ref.f[0], hst, hsr, hstd,
                         hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta,
                         pathtype);

    const double hts = profile.h.front() + ref.htg;
    const double hrs = profile.h.back() + ref.hrg;

    const double omega = path_fraction_sea(d, zone, 3);

    EXPECT_NEAR(ae, ref.ae_ref, kTol) << "ae";
    EXPECT_NEAR(dtot, ref.dtot_ref, kTol) << "dtot";
    EXPECT_NEAR(hts, ref.hts_ref, kTol) << "hts";
    EXPECT_NEAR(hrs, ref.hrs_ref, kTol) << "hrs";
    EXPECT_NEAR(theta_t, ref.theta_t_ref, kTol) << "theta_t";
    EXPECT_NEAR(theta_r, ref.theta_r_ref, kTol) << "theta_r";
    EXPECT_NEAR(theta, ref.theta_ref, kTol) << "theta";
    EXPECT_NEAR(hm, ref.hm_ref, kTol) << "hm";
    EXPECT_NEAR(hte, ref.hte_ref, kTol) << "hte";
    EXPECT_NEAR(hre, ref.hre_ref, kTol) << "hre";
    EXPECT_NEAR(hstd, ref.hstd_ref, kTol) << "hstd";
    EXPECT_NEAR(hsrd, ref.hsrd_ref, kTol) << "hsrd";
    EXPECT_NEAR(dlt, ref.dlt_ref, kTol) << "dlt";
    EXPECT_NEAR(dlr, ref.dlr_ref, kTol) << "dlr";
    EXPECT_EQ(pathtype, ref.pathtype_ref) << "pathtype";
    EXPECT_NEAR(dtm, ref.dtm_ref, kTol) << "dtm";
    EXPECT_NEAR(dlm, ref.dlm_ref, kTol) << "dlm";
    EXPECT_NEAR(b0, ref.b0_ref, kTol) << "b0";
    EXPECT_NEAR(omega, ref.omega_ref, kTol) << "omega";
    EXPECT_NEAR(DN, ref.DN_ref, kTol) << "DN";
    EXPECT_NEAR(N0, ref.N0_ref, kTol) << "N0";

    // --- Per (f, p) transmission-loss terms ---

    const double d3D = std::sqrt(dtot * dtot + std::pow((hts - hrs) / 1000.0, 2.0));

    for (size_t i = 0; i < ref.f.size(); ++i) {
      SCOPED_TRACE("Row " + std::to_string(i) + " (f=" + std::to_string(ref.f[i]) +
                   " GHz, p=" + std::to_string(ref.p[i]) + "%)");

      const double f = ref.f[i];
      const double p = ref.p[i];

      double Lbfsg, Lb0p, Lb0b;
      pl_los(d3D, f, p, b0, omega, ref.temp, ref.press, dlt, dlr, Lbfsg, Lb0p,
            Lb0b);

      const double Lba = tl_anomalous(dtot, dlt, dlr, ref.dct, ref.dcr, dlm,
                                      hts, hrs, hte, hre, hm, theta_t, theta_r,
                                      f, p, ref.temp, ref.press, omega, ae, b0);

      // Bullington diffraction loss for the actual profile ...
      const double Lbulla = dl_bull(d, g, hts, hrs, ae, f);

      // ... and again for the equivalent smooth path, to get the
      // spherical-Earth diffraction loss (eq. 39/40 of P.452).
      const double hts1 = hts - hstd;
      const double hrs1 = hrs - hsrd;
      std::vector<double> h1(n, 0.0);
      const double Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);
      (void)Lbulla;
      (void)Lbulls;  // Not compared directly - only Ldsph (below) is.

      const double hte1 = hts1;
      const double hre1 = hrs1;
      const double Ldsph = dl_se(dtot, hte1, hre1, ae, f, omega, ref.pol);

      double Ldp, Ld50;
      dl_p(d, g, hts, hrs, hstd, hsrd, f, omega, p, b0, DN, ref.pol, Ldp, Ld50);

      const double Lbs = tl_tropo(dtot, theta, f, p, ref.temp, ref.press, N0,
                                  ref.Gt, ref.Gr);

      // tl_p452 mutates its d/h/g/zone arguments in place (it re-applies
      // the Step 4 clutter override), so pass fresh copies each call.
      std::vector<double> dCall = profile.d;
      std::vector<double> hCall = profile.h;
      std::vector<double> gCall = g;
      std::vector<int> zoneCall = profile.zone;
      const double Lb =
          tl_p452(maps, f, p, dCall, hCall, gCall, zoneCall, ref.htg, ref.hrg,
                 ref.phit_e, ref.phit_n, ref.phir_e, ref.phir_n, ref.Gt, ref.Gr,
                 ref.pol, ref.dct, ref.dcr, ref.press, ref.temp);

      EXPECT_NEAR(Lbfsg, ref.Lbfsg_ref[i], kTol) << "Lbfsg";
      EXPECT_NEAR(Lb0p, ref.Lb0p_ref[i], kTol) << "Lb0p";
      EXPECT_NEAR(Lb0b, ref.Lb0b_ref[i], kTol) << "Lb0b";
      EXPECT_NEAR(Ldsph, ref.Ldsph_ref[i], kTol) << "Ldsph";
      EXPECT_NEAR(Ld50, ref.Ld50_ref[i], kTol) << "Ld50";
      EXPECT_NEAR(Ldp, ref.Ldp_ref[i], kTol) << "Ldp";
      EXPECT_NEAR(Lbs, ref.Lbs_ref[i], kTol) << "Lbs";
      EXPECT_NEAR(Lba, ref.Lba_ref[i], kTol) << "Lba";
      EXPECT_NEAR(Lb, ref.Lb_ref[i], kTol) << "Lb";
    }
  }
}
