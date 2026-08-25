#!/usr/bin/env bash

set -euo pipefail

# Determine the TALYS installation directory, independently of where
# the script is called from.

talys_dir=$(cd "$(dirname "$0")" && pwd)
source_dir="$talys_dir/source"

# Verify that the expected directories and build files exist.

if [[ ! -d "$source_dir" ]]; then
  echo "TALYS installation error: source directory not found:" >&2
  echo "  $source_dir" >&2
  exit 1
fi

if [[ ! -f "$source_dir/Makefile" ]]; then
  echo "TALYS installation error: Makefile not found:" >&2
  echo "  $source_dir/Makefile" >&2
  exit 1
fi

if [[ ! -d "$talys_dir/structure" ]]; then
  echo "TALYS installation error: structure database not found:" >&2
  echo "  $talys_dir/structure" >&2
  exit 1
fi

echo
echo "Installing TALYS"
echo "Installation directory: $talys_dir"
echo

# Pass all command-line arguments directly to make. This permits, e.g.:
#
# ./install_talys.bash FC=ifx
# ./install_talys.bash FFLAGS="-O3 -march=native"
# ./install_talys.bash FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"   (the optimal choice for MacOS)
# ./install_talys.bash clean

make -C "$source_dir" clean
make -C "$source_dir" "$@"

talys_exe="$talys_dir/bin/talys"

if [[ -x "$talys_exe" ]]; then
  echo
  echo "TALYS executable:"
  echo "  $talys_exe"
  echo
  echo "If not already done, add the following lines to your shell configuration:"
  echo
  echo "  export TALYS_DIR=\"$talys_dir\""
  echo "  export PATH=\"\$TALYS_DIR/bin:\$PATH\""
  echo "  export TALYS_USER=\"Your Name\""
  echo
  echo "Alternatively, edit code_dir in source/machine.f90 and rebuild TALYS."
  echo
fi
