#! /bin/bash

ncpus=1
if test -x /usr/bin/getconf; then
    ncpus=$(/usr/bin/getconf _NPROCESSORS_ONLN)
fi

echo "Timestamp" && date
mkdir cmake-build
cd cmake-build

OSX_SDK="/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk"
#      -DNETGEN_SUPPORT=ON \

#      OCE_INCLUDE_PATH
#      OCE_LIB_PATH

#      -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE \


#      -DCMAKE_MACOSX_RPATH=OFF \
#      -DCMAKE_SKIP_RPATH=ON \
#      -DBUILD_WITH_INSTALL_RPATH=OFF \
#      -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=On

#      -DCMAKE_CXX_COMPILER=$PREFIX/bin/c++ \

cmake -DCMAKE_BUILD_TYPE:STRING=Release \
      -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
      -DCMAKE_C_COMPILER=/usr/bin/clang \
      -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DCMAKE_OSX_SYSROOT=$OSX_SDK \
      -DCMAKE_OSX_DEPLOYMENT_TARGET=10.9 \
      -DBoost_DIR=$PREFIX/share/cmake-3.0/Modules/ \
      -DBoost_INCLUDE_DIR=$PREFIX/include \
      -DOCE_DIR=$PREFIX/lib \
      -DCMAKE_OSX_DEPLOYMENT_TARGET= \
      -DCMAKE_OSX_SYSROOT= \
      -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DSMESH_TESTING=OFF \
      -DCMAKE_MACOSX_RPATH=OFF \
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