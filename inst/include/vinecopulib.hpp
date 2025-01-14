// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
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

// silences all the BOOST_CONCEPT warnings bullshit
// https://stackoverflow.com/questions/13930894/how-to-disable-boost-concept-check
#include <boost/concept/assert.hpp>
#undef BOOST_CONCEPT_ASSERT
#define BOOST_CONCEPT_ASSERT(Model)
#include <boost/concept_check.hpp>

#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/vinecop/class.hpp>
#include <wdm/eigen.hpp>
