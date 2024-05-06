#!/usr/bin/env bash

set -e
set -u

# arguments
[ "${1:-}" = "clean" ] && clean=true || clean=false

# verify pythia8 is installed
which pythia8-config || (echo "ERROR: pythia8 is not found" >&2 && exit 1)

# cd to the stringspinner dir
src_dir=$(cd $(dirname ${BASH_SOURCE[0]:-$0})/stringspinner && pwd -P)
pushd $src_dir

# apply patches
for patch in ../patches/*.patch; do
  echo "[+] apply patch '$(basename $patch)'"
  git apply $patch
done

# build
echo "[+] building stringspinner"
./configure $(pythia8-config --prefix)
if $clean; then make clean; fi
make

# revert patches, since they are only needed for building
for patch in ../patches/*.patch; do
  echo "[+] revert patch '$(basename $patch)'"
  git apply -R --whitespace=nowarn $patch
done

popd
echo "[+] stringspinner built"
