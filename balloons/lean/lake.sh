#!/usr/bin/env bash
set -euo pipefail
balloons_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
balloons_bin="$balloons_root/.toolchain/lean-4.35.0-linux/bin"
export PATH="$balloons_bin:$PATH"
export MATHLIB_CACHE_DIR="$balloons_root/.lake/cache"
cd "$balloons_root"
exec "$balloons_bin/lake" "$@"
