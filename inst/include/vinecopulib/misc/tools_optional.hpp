// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#if __cplusplus >= 201703L
#include <optional>
namespace vinecopulib {
namespace optional {

template<typename T>
using optional = std::optional<T>;

template<typename T>
bool
has_value(const optional<T>& opt)
{
  return opt.has_value();
}

template<typename T>
T&
value(optional<T>& opt)
{
  return opt.value();
}

template<typename T>
const T&
value(const optional<T>& opt)
{
  return opt.value();
}

} // namespace optional
} // namespace vinecopulib
#else
#include <boost/optional.hpp>
namespace vinecopulib {
namespace optional {

template<typename T>
using optional = boost::optional<T>;

template<typename T>
bool
has_value(const optional<T>& opt)
{
  return static_cast<bool>(opt);
}

template<typename T>
T&
value(optional<T>& opt)
{
  return opt.get();
}

template<typename T>
const T&
value(const optional<T>& opt)
{
  return opt.get();
}

} // namespace optional
} // namespace vinecopulib
#endif