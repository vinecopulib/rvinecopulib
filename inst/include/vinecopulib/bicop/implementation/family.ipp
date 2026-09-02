// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>
#include <unordered_map>

namespace vinecopulib {

namespace {

// The two directions are kept as separate maps rather than a bidirectional
// container: the table is fixed and tiny, and both lookups are on hot paths
// (family names appear in every serialized model).
inline const std::vector<std::pair<BicopFamily, std::string>>&
family_name_table()
{
  static const std::vector<std::pair<BicopFamily, std::string>> table = {
    { BicopFamily::indep, "Independence" },
    { BicopFamily::gaussian, "Gaussian" },
    { BicopFamily::student, "Student" },
    { BicopFamily::clayton, "Clayton" },
    { BicopFamily::gumbel, "Gumbel" },
    { BicopFamily::frank, "Frank" },
    { BicopFamily::joe, "Joe" },
    { BicopFamily::bb1, "BB1" },
    { BicopFamily::bb6, "BB6" },
    { BicopFamily::bb7, "BB7" },
    { BicopFamily::bb8, "BB8" },
    { BicopFamily::tawn, "Tawn" },
    { BicopFamily::tll, "TLL" }
  };
  return table;
}

struct BicopFamilyHash
{
  size_t operator()(BicopFamily family) const noexcept
  {
    return std::hash<int>()(static_cast<int>(family));
  }
};

inline const std::unordered_map<BicopFamily, std::string, BicopFamilyHash>&
enum_to_name()
{
  static const std::unordered_map<BicopFamily, std::string, BicopFamilyHash>
    map(family_name_table().begin(), family_name_table().end());
  return map;
}

inline const std::unordered_map<std::string, BicopFamily>&
name_to_enum()
{
  static const std::unordered_map<std::string, BicopFamily> map = [] {
    std::unordered_map<std::string, BicopFamily> m;
    for (const auto& e : family_name_table()) {
      m.emplace(e.second, e.first);
    }
    return m;
  }();
  return map;
}
}

//! @brief Converts a BicopFamily into a string with its name.
//! @param family The family.
inline std::string
get_family_name(BicopFamily family)
{
  return enum_to_name().at(family);
}

//! @brief Converts a string name into a BicopFamily.
//! @param family The family name.
inline BicopFamily
get_family_enum(const std::string& family)
{
  return name_to_enum().at(family);
}
}
