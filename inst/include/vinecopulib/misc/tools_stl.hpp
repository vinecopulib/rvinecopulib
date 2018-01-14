// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cstddef>
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <limits>

namespace vinecopulib {

namespace tools_stl {

using namespace std;

template<class... T>
void unused(T &&...)
{

}

template<typename T, typename V>
V unaryFunc_or_nan(T f, V y)
{
    if ((boost::math::isnan)(y)) {
        return std::numeric_limits<V>::quiet_NaN();
    } else {
        return f(y);
    }
}

template<typename T, typename V>
V binaryFunc_or_nan(T f, V u1, V u2)
{
    if ((boost::math::isnan)(u1) | (boost::math::isnan)(u2)) {
        return std::numeric_limits<V>::quiet_NaN();
    } else {
        return f(u1, u2);
    }
}

template<typename T>
std::vector <size_t> get_order(const std::vector <T> &x)
{
    std::vector <size_t> order(x.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(
        order.begin(),
        order.end(),
        [&](size_t i, size_t j) -> bool { return (x[i] < x[j]); }
    );
    return order;
}

template<typename T>
bool is_member(T element, vector <T> set)
{
    return find(set.begin(), set.end(), element) != set.end();
}

template<class T>
vector <T> intersect(vector <T> x, vector <T> y)
{
    sort(x.begin(), x.end());
    sort(y.begin(), y.end());
    vector <T> common;
    set_intersection(
        x.begin(), x.end(),
        y.begin(), y.end(),
        back_inserter(common)
    );

    return common;
}

template<class T>
size_t find_position(T x, vector <T> vec)
{
    return distance(vec.begin(), find(vec.begin(), vec.end(), x));
}

template<class T>
vector <T> set_diff(vector <T> x, vector <T> y)
{
    sort(x.begin(), x.end());
    sort(y.begin(), y.end());
    vector <T> different;
    set_difference(
        x.begin(), x.end(),
        y.begin(), y.end(),
        back_inserter(different)
    );

    return different;
}

template<class T>
vector <T> cat(vector <T> x, const vector <T> &y)
{
    x.reserve(x.size() + y.size());
    x.insert(x.end(), y.begin(), y.end());
    return x;
}

template<class T>
vector <T> cat(T x, const vector <T> &y)
{
    vector <T> out(1);
    out[0] = x;
    out.reserve(1 + y.size());
    out.insert(out.end(), y.begin(), y.end());
    return out;
}

template<class T>
vector <T> set_sym_diff(vector <T> x, vector <T> y)
{
    vector <T> dxy = set_diff(x, y);
    auto dyx = set_diff(y, x);
    return cat(dxy, dyx);
}

template<class T>
void reverse(vector <T> &x)
{
    reverse(x.begin(), x.end());
}

template<class T>
bool is_same_set(vector <T> x, vector <T> y)
{
    auto z = intersect(x, y);
    return ((z.size() == x.size()) & (z.size() == y.size()));
}

//! Integer sequence starting at 1
inline vector <size_t> seq_int(size_t from, size_t length)
{
    vector <size_t> seq(length);
    iota(seq.begin(), seq.end(), from);
    return seq;
}
}

}
