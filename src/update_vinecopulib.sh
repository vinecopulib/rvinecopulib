#!/usr/bin/env bash

set -euo pipefail

backend_ref="${1:-main}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
include_dir="$script_dir/../inst/include"
work_dir="$(mktemp -d)"
backend_dir="$work_dir/vinecopulib"
backend_url="git@github.com:vinecopulib/vinecopulib.git"

cleanup() {
  rm -rf -- "$work_dir"
}
trap cleanup EXIT

git init --quiet "$backend_dir"
git -C "$backend_dir" remote add origin "$backend_url"
git -C "$backend_dir" fetch --quiet --depth 1 origin "$backend_ref"
git -C "$backend_dir" checkout --quiet --detach FETCH_HEAD

backend_commit="$(git -C "$backend_dir" rev-parse HEAD)"

# Keep package-owned headers such as vinecopulib-wrappers.hpp, but replace the
# complete upstream header tree, including documentation headers.
rm -rf -- "$include_dir/vinecopulib"
cp -R "$backend_dir/include/." "$include_dir/"
cp "$backend_dir/LICENSE" "$include_dir/vinecopulib/LICENSE"
printf '%s\n' "$backend_commit" > "$include_dir/vinecopulib/REVISION"
find "$include_dir/vinecopulib" -type f -exec chmod 0644 {} +

printf 'Imported vinecopulib %s at %s\n' "$backend_ref" "$backend_commit"
