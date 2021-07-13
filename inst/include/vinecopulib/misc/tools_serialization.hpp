// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <vinecopulib/misc/nlohmann_json.hpp>
#include <vinecopulib/misc/triangular_array.hpp>

namespace vinecopulib {

namespace tools_serialization {

//! conversion from Eigen::Matrix to nlohmann::json
//!
//! @param matrix The Eigen::Matrix to convert.
//! @return the corresponding nlohmann::json.
template<class T>
inline nlohmann::json
matrix_to_json(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix)
{

  nlohmann::json output;
  output["shape"] = { matrix.rows(), matrix.cols() };
  nlohmann::json json_data;
  auto it = matrix.data();
  auto last = matrix.data() + matrix.size();
  for (; it != last; ++it) {
    json_data.push_back(*it);
  }
  output["data"] = json_data;

  return output;
}

//! conversion from vinecopulib::TriangularArray to nlohmann::json
//!
//! @param array The vinecopulib::TriangularArray to convert.
//! @return the corresponding nlohmann::json.
template<class T>
inline nlohmann::json
triangular_array_to_json(const TriangularArray<T>& array)
{
  nlohmann::json output;
  size_t d = array.get_dim();
  size_t trunc_lvl = array.get_trunc_lvl();
  output["d"] = d;
  output["t"] = trunc_lvl;

  nlohmann::json json_data;
  for (size_t i = 0; i < std::min(d - -1, trunc_lvl); i++) {
    nlohmann::json row;
    for (size_t j = 0; j < d - 1 - i; j++) {
      row.push_back(array(i, j));
    }
    json_data.push_back(row);
  }
  output["data"] = json_data;

  return output;
}

//! conversion from std::vector to nlohmann::json
//!
//! @param vec The std::vector to convert.
//! @return the corresponding nlohmann::json.
template<class T>
inline nlohmann::json
vector_to_json(const std::vector<T>& vec)
{
  nlohmann::json output = vec;
  return output;
}

//! conversion from nlohmann::json to Eigen::Matrix
//!
//! @param input The nlohmann::json to convert.
//! @return the corresponding Eigen::Matrix.
template<typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
json_to_matrix(const nlohmann::json& input)
{

  size_t rows = input["shape"][0];
  size_t cols = input["shape"][1];
  std::vector<double> vec = input["data"];

  Eigen::MatrixXd matrix;
  matrix = Eigen::MatrixXd::Map(&vec[0], rows, cols);

  return matrix.cast<T>();
}

//! conversion from nlohmann::json to vinecopulib::TriangularArray
//!
//! @param input The nlohmann::json to convert.
//! @return the corresponding vinecopulib::TriangularArray
template<typename T>
inline TriangularArray<T>
json_to_triangular_array(const nlohmann::json& input)
{

  std::vector<std::vector<T>> vec = input["data"];
  return TriangularArray<T>(vec);
}

//! conversion from nlohmann::json to std::vector
//!
//! @param input The nlohmann::json to convert.
//! @return the corresponding std::vector.
template<typename T>
inline std::vector<T>
json_to_vector(const nlohmann::json& input)
{

  std::vector<T> res = input;
  return res;
}

inline nlohmann::json
file_to_json(const std::string& filename)
{
  nlohmann::json output;
  std::ifstream file(filename);
  file >> output;
  return output;
}

inline void
json_to_file(const std::string& filename, const nlohmann::json& json)
{
  std::ofstream file(filename);
  file << json << std::endl;
}
}
}
