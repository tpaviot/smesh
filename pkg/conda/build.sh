#! /bin/bash

set -e

ncpus=4

echo "Timestamp" && date
mkdir cmake-build
cd cmake-build
cmake -DCMAKE_BUILD_TYPE:STRING=Release \
      -DOCE_DIR=$PREFIX/lib \
      -DCMAKE_INSTALL_PREFIX=$PREFIX \
      ..
echo ""
echo "Timestamp" && date
# make for py3
echo "Starting build with -j$ncpus ..."
# travis-ci truncates when there are more than 10,000 lines of output.
# trim them to see test results.
make -j$ncpus install | grep Built
# Run tests
echo "Timestamp" && date