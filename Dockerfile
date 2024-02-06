#------------------------------------------------------------------------------
# CSIRO EMSr Image Build Script
#------------------------------------------------------------------------------
# The recommended base image is one of the onaci/ereefs-netcdf-base variants,
# as those have been designed to include all the EMS dependencies.
# Allow it to be overridden in order to choose *which* variant, or even
# something completely different.
ARG BASE_IMAGE="onaci/ereefs-netcdf-base:ems"
FROM ${BASE_IMAGE}

# Record the actual base image used from the FROM command as a label.
ARG BASE_IMAGE
LABEL org.opencontainers.image.base.name=${BASE_IMAGE}

# Enable Bash in RUN commands, and ensure that any commands with
# pipes exit on the first failure.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Prepare a base directory for all the EMS code
ARG EMS_BASE="/usr/local/src/ems"
ENV EMS_BASE="${EMS_BASE}"
WORKDIR ${EMS_BASE}

# Install all the EMS Source Code
COPY ./ ./

# Selectively remove some subdirectories that we want to omit
# (List can be overridden by a build argument)
ARG EMS_RM_DIRS=""
RUN if [ -n "${EMS_RM_DIRS}" ]; then for d in $EMS_RM_DIRS; do rm -rf "${EMS_BASE}/${d}"; done; fi

# Compile the EMS components
# We build with OpenMP support by default for use in a docker container
# Different options can be specified by overriding the EMS_CONFIGURE_OPTS build argument
# (The base image *does* include NetCDF and HDF5 with Open MPI support, so that is an option)
ENV EMS_BUILD_LOG="${EMS_BASE}/ems_build.log"
ARG EMS_CONFIGURE_OPTS="--enable-omp"
RUN make distclean || true \
    && conf/configure ${EMS_CONFIGURE_OPTS} 2>&1 | tee -a "${EMS_BUILD_LOG}"
RUN make clean \
    && make 2>&1 | tee -a "${EMS_BUILD_LOG}"
RUN make check install 2>&1 | tee -a "${EMS_BUILD_LOG}"

# Symlink the EMS executable(s) into the default path
# (This will fail if the executable did not build for some reason)
ENV SHOC="${EMS_BASE}/model/hd/shoc" \
    COMPAS="${EMS_BASE}/model/hd-us/compas"
RUN if [ -n "${SHOC}" ]; then chmod 0755 "${SHOC}" && ln -s "${SHOC}" /usr/local/bin/shoc; fi; \
    if [ -n "${COMPAS}" ]; then chmod 0755 "${COMPAS}" && ln -s "${COMPAS}" /usr/local/bin/compas; fi; \
    ln -s "${EMS_BASE}/model/tests/run-tests.sh" /usr/local/bin/run-ems-tests

# Optionally run all available unit tests
# (Override the EMS_TEST_RUN build argument to have a value of 0 to skip testing at build time)
ARG EMS_TEST_RUN=1
ARG EMS_TEST_INCLUDE="${EMS_BASE}/model/tests"
ENV EMS_TEST_INCLUDE="${EMS_TEST_INCLUDE}"
ARG EMS_TEST_EXCLUDE=""
ENV EMS_TEST_EXCLUDE="${EMS_TEST_EXCLUDE}"
ENV EMS_TEST_LOG="${EMS_BASE}/model/tests/ems_test.log"
RUN if [ $EMS_TEST_RUN -eq 1 ]; then run-ems-tests; fi

# Encode EMS metadata in labels
LABEL au.csiro.ems.base=${EMS_BASE} \
      au.csiro.ems.shoc=${SHOC} \
      au.csiro.ems.compas=${COMPAS}

# Configure the default entrypoint to be a default EMS executable
ENTRYPOINT ["/bin/bash", "-c" ]
CMD [ "${SHOC:-$COMPAS}", "-v"]
