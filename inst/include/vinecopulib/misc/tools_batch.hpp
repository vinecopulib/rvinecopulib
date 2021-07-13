// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <cmath>
#include <vector>

namespace vinecopulib {

namespace tools_batch {

struct Batch
{
  size_t begin;
  size_t size;
};

inline size_t
compute_num_batches(size_t num_tasks, size_t num_threads)
{
  if (num_tasks < num_threads)
    return num_tasks;
  size_t num_batches =
    num_threads *
    static_cast<size_t>(std::floor(std::sqrt(num_tasks / num_threads)));
  return std::min(num_tasks, num_batches);
}

inline std::vector<Batch>
create_batches(size_t num_tasks, size_t num_threads)
{
  if (num_tasks == 0)
    return { Batch{ 0, 0 } };
  num_threads = std::max(static_cast<size_t>(1), num_threads);

  size_t num_batches = compute_num_batches(num_tasks, num_threads);
  std::vector<Batch> batches(num_batches);

  size_t min_size = num_tasks / num_batches;
  ptrdiff_t rem_size = num_tasks % num_batches;
  for (size_t i = 0, k = 0; i < num_tasks; k++) {
    batches[k] = Batch{ i, min_size + (rem_size-- > 0) };
    i += batches[k].size;
  }

  return batches;
}
}
}
