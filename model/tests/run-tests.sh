#!/bin/bash
#
# Runs the EMS test suite(s)
#
set +e

SUCCESS_LIST=""
SUCCESS_COUNT=0

FAIL_LIST=""
FAIL_COUNT=0

echo "SHOC = '${SHOC}' => '$(which shoc)'"
echo "COMPAS = '${COMPAS}' => '$(which compas)'"
echo
echo "CURL_VERSION = '${CURL_VERSION}'"
echo "DAP_VERSION = '${DAP_VERSION}'"
echo "GDAL_VERSION = '${GDAL_VERSION}'"
echo "HDF5_VERSION = '${HDF5_VERSION}'"
echo "NETCDF_VERSION = '${NETCDF_VERSION}'"
echo "NCO_VERSION = '${NCO_VERSION}'"
echo "PROJ_VERSION = '${PROJ_VERSION}'"
echo

EMS_TEST_INCLUDE="${EMS_TEST_INCLUDE:-${EMS_BASE}/model/tests}"
echo "EMS_TEST_INCLUDE = '${EMS_TEST_INCLUDE}'"

EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE:-} ${EMS_BASE}/model/tests/run-tests.sh"
if [ -z "${SHOC}" ]; then
    EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE} '${EMS_BASE}/model/tests/hd/*'"
else
    EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE} '${EMS_BASE}/model/tests/hd/run_tests'"    
fi
if [ -z "${COMPAS}" ]; then
    EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE} '${EMS_BASE}/model/tests/hd-us/*'"
else
    EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE} '${EMS_BASE}/model/tests/hd-us/run_tests'"
fi
echo "EMS_TEST_EXCLUDE = '${EMS_TEST_EXCLUDE}'"

TEST_COUNT=0

for ti in $EMS_TEST_INCLUDE; do
    FIND_CMD="find ${ti} -type f \( -regex '^.*/run_test[0-9]*$' -o -name 'run_all.sh' -o -regex '^.*/run_[^\.]*$' \) ! -iname '.*'"
    for te in $EMS_TEST_EXCLUDE; do
        FIND_CMD="${FIND_CMD} ! -path ${te}"
    done
    echo "FIND_CMD = '${FIND_CMD}'"
    TEST_FILES=$(eval "${FIND_CMD}")
    if [ -z "${TEST_FILES}" ]; then
        echo "There are no tests to run: Test FAILED"
        FAIL_LIST="${FAIL_LIST} ${ti}"
        FAIL_COUNT=$((FAIL_COUNT+1))
    else
        echo "Test files to run:"
        echo "${TEST_FILES}"
        for tf in $TEST_FILES; do
            chmod 0755 "${tf}"
            cd $(dirname "${tf}")
            echo
            echo "#-----------------------------------------------------"
            echo "# Running test file '${tf}' from $(pwd)..."
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
exit $FAIL_COUNT
