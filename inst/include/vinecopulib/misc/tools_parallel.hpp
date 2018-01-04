// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_interface.hpp>

namespace tools_parallel {

//! maps a function on a list of items, possibly running tasks in parallel.
//! @param f function to be mapped.
//! @param items an objects containing the items on which `f` shall be mapped;
//!     must allow for `auto` loops (i.e., `std::begin(I)`/`std::end(I)` must be 
//!     defined).
//! @param num_threads the number of concurrent threads to use in the pool.
template<class F, class I>
void map_on_pool(F &&f, I &&items, size_t num_threads)
{
    if (num_threads <= 1) {
        for (const auto &item : items) {
            f(item);
        }
    } else {
        ThreadPool pool(num_threads);
        for (const auto &item : items) {
            pool.push(f, item);
        }
        pool.join();
    }
}

}
