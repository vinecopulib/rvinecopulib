#!/usr/bin/env bash

set -euo pipefail

backend_ref="${1:-main}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
include_dir="$script_dir/../inst/include"
work_dir="$(mktemp -d)"
backend_dir="$work_dir/vinecopulib"

cleanup() {
  rm -rf -- "$work_dir"
}
trap cleanup EXIT

git clone \
  --depth 1 \
  --branch "$backend_ref" \
  --single-branch \
  git@github.com:vinecopulib/vinecopulib.git \
  "$backend_dir"

backend_commit="$(git -C "$backend_dir" rev-parse HEAD)"

# Keep package-owned headers such as vinecopulib-wrappers.hpp, but replace the
# complete upstream header tree, including documentation headers.
rm -rf -- "$include_dir/vinecopulib"
cp -R "$backend_dir/include/." "$include_dir/"
cp "$backend_dir/LICENSE" "$include_dir/vinecopulib/LICENSE"

printf 'Imported vinecopulib %s at %s\n' "$backend_ref" "$backend_commit"
