#!/bin/bash
#
# Runs the EMS test suite(s)
#

# Don't exit for errors: we will handle (and count) those.
set +e

# Send all output to a logfile as well as stdout / stderr
EMS_TEST_LOG="${EMS_TEST_LOG:-ems_test.log}"
if [ -f "${EMS_TEST_LOG}" ]; then
    truncate -s 0 "${EMS_TEST_LOG}"
fi
exec > >(tee -a "${EMS_TEST_LOG}")
exec 2> >(tee -a "${EMS_TEST_LOG}" >&2)

# Dump some useful information about the test environment
# (note: some of these vars are set by the base image)
echo "Testing EMS at $(date --rfc-3339=seconds)"
echo
echo "SHOC = '${SHOC}' => '$(which shoc)'"
echo "COMPAS = '${COMPAS}' => '$(which compas)'"
echo
echo "CURL_VERSION = '${CURL_VERSION:-}', curl binary = '$(which curl)' => $(curl --version | head -n 1)"
echo "DAP_VERSION = '${DAP_VERSION:-}', getdap4 binary = '$(which getdap4)' => $(getdap4 -V 2>&1)"
echo "GDAL_VERSION = '${GDAL_VERSION:-}', gdalinfo binary = '$(which gdalinfo)', => $(gdalinfo --version 2>&1)"
echo "HDF5_VERSION = '${HDF5_VERSION:-}', h5dump binary = '$(which h5dump)' => $(h5dump --version)"
echo "NETCDF_VERSION = '${NETCDF_VERSION:-}', ncdump binary = '$(which ncdump)' => $(ncdump 2>&1 | tail -n 1)"
echo "NCO_VERSION = '${NCO_VERSION:-}', ncks binary = '$(which ncks)' => $(ncks --version 2>&1 | xargs)"
echo "PROJ_VERSION = '${PROJ_VERSION:-}', proj binary = '$(which proj)' => $(proj 2>&1 | head -n 1)"
echo

# Initialise vars to hold test result statistics
FAIL_LIST=""
FAIL_COUNT=0
SUCCESS_LIST=""
SUCCESS_COUNT=0
TEST_COUNT=0

# Parse the include and exclude variables used to select target tests
EMS_TEST_INCLUDE="${EMS_TEST_INCLUDE:-${EMS_BASE}/model/tests}"
echo "EMS_TEST_INCLUDE = '${EMS_TEST_INCLUDE}'"

EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE:-} ${EMS_BASE}/model/tests/run-tests.sh '${EMS_BASE}/model/tests/hd/run_tests'"
if [ ! -f "${COMPAS}" ]; then
    # Valid scenario: skip the COMPAS tests
    EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE} '${EMS_BASE}/model/tests/hd-us/*'"
else
    EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE} '${EMS_BASE}/model/tests/hd-us/run_tests'"
fi
echo "EMS_TEST_EXCLUDE = '${EMS_TEST_EXCLUDE}'"


for ti in $EMS_TEST_INCLUDE; do
    FIND_CMD="find ${ti} -type f \( -regex '^.*/run_test[0-9]*$' -o -name 'run_all.sh' -o -regex '^.*/run_[^\.]*$' \) ! -iname '.*'"
    for te in $EMS_TEST_EXCLUDE; do
        FIND_CMD="${FIND_CMD} ! -path ${te}"
    done
    echo "FIND_CMD = '${FIND_CMD}'"
    TEST_FILES=$(eval "${FIND_CMD}")
    if [ -z "${TEST_FILES}" ]; then
        echo "There are no test files to run below include path '${ti}'"
    else
        echo "Test files to run below include path '${ti}':"
        echo "${TEST_FILES}"
        for tf in $TEST_FILES; do
            chmod 0755 "${tf}"
            cd $(dirname "${tf}")
            echo
            echo "#-----------------------------------------------------"
            echo "# Running test file '${tf}'"
            echo "# CWD = '$(pwd)'"
            echo "#-----------------------------------------------------"
            ${tf}
            TEST_RESULT=$?
            echo "#......................................................"
            if [ $TEST_RESULT -eq 0 ]; then
                echo "exitcode: ${TEST_RESULT} => Test SUCCEEDED"
                SUCCESS_LIST="${SUCCESS_LIST} ${tf}"
                SUCCESS_COUNT=$((SUCCESS_COUNT+1))
            else
                echo "exitcode: ${TEST_RESULT} => Test FAILED"
                FAIL_LIST="${FAIL_LIST} ${tf}"
                FAIL_COUNT=$((FAIL_COUNT+1))
            fi
            TEST_COUNT=$((TEST_COUNT+1))
        done
    fi
done

echo
echo "#-----------------------------------------------------"
echo
echo "${SUCCESS_COUNT} of ${TEST_COUNT} tests passed:"
echo "${SUCCESS_LIST}" | tr ' ' '\n'
echo
echo "${FAIL_COUNT} of ${TEST_COUNT} tests failed:"
echo "${FAIL_LIST}" | tr ' ' '\n'
if [ $TEST_COUNT -eq 0 ]; then
    exit -1
else
    exit $FAIL_COUNT
fi
