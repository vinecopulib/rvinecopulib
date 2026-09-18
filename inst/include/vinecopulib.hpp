// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#ifndef BOOST_NO_AUTO_PTR
#define BOOST_NO_AUTO_PTR
#endif

#ifndef BOOST_ALLOW_DEPRECATED_HEADERS
#define BOOST_ALLOW_DEPRECATED_HEADERS
#endif

#ifndef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#else
#undef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#endif

#ifndef BOOST_ALL_NO_LIB
#define BOOST_ALL_NO_LIB
#endif

// wdm/random.hpp selects Boost's mt19937 over std::mt19937 on this; keep it in
// sync with VINECOPULIB_DEFINITIONS so header-only and precompiled agree.
#ifndef USE_BOOST
#define USE_BOOST
#endif

#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/vinecop/class.hpp>
#include <wdm/eigen.hpp>
