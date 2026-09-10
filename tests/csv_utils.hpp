#pragma once

// Minimal CSV reader mirroring matlab/private/readcsv.m: reads a
// comma-separated file and returns the data rows as strings, skipping the
// header row. Numeric conversion is left to the caller (via std::stod /
// std::stoi), just like the MATLAB script uses str2double on demand.

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

inline std::vector<std::vector<std::string>> readCsv(const std::string& path) {
  std::ifstream in(path);
  if (!in.is_open()) {
    throw std::runtime_error("File " + path + " cannot be found.");
  }

  std::vector<std::vector<std::string>> rows;
  std::string line;
  bool isHeader = true;

  while (std::getline(in, line)) {
    if (isHeader) {
      // Skip the header row, exactly as readcsv.m does.
      isHeader = false;
      continue;
    }
    if (line.empty()) {
      continue;
    }

    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) {
      fields.push_back(field);
    }
    rows.push_back(std::move(fields));
  }

  return rows;
}
