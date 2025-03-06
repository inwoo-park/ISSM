#!/bin/bash
set -eu


## Constants
#
VER="3.20.5"

PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src" # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install" # Set to location where external package should be installed

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/petsc-lite-${VER}.tar.gz" "petsc-${VER}.tar.gz"

# Unpack source
tar -zxvf petsc-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} ${PETSC_DIR}
mkdir -p ${PETSC_DIR}

# Move source to $PETSC_DIR
mv petsc-${VER}/* ${PETSC_DIR}
rm -rf petsc-${VER}

# Configure
#
# NOTE:
# - Added -Wno-error=implicit-function-declaration to CFLAGS for Clang >= 12. 
#	(may need to remove it for earlier versions not using the C99 standard).
#
cd ${PETSC_DIR}
./configure \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-pic=1 \
	--download-fblaslapack=1 \
	--download-metis=1 \
	--download-mpich="https://www.mpich.org/static/downloads/4.2.0/mpich-4.2.0.tar.gz" \
	--download-mumps=1 \
	--download-parmetis=1 \
	--download-scalapack=1 \
	--download-zlib=1

# Compile and install
make
make install
