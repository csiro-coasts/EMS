#!/bin/bash -e

THIS_SCRIPT="${BASH_SOURCE[0]}"
THIS_DIR=$( cd "$(dirname "${THIS_SCRIPT}")" && pwd)

"${THIS_DIR}/run.sh" "test_basic"

"${THIS_DIR}/run.sh" "test_advanced" "sedi"
