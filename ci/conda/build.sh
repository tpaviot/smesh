# Configure step
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=$PREFIX \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_PREFIX_PATH=$PREFIX \
 -DCMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
 -DBoost_INCLUDE_DIRS=$PREFIX/include/ \
 -DSMESH_TESTING=ON \
 -DSMESH_BUILD_NETGENPLUGIN=ON \
 -DNETGEN_NGLIB_INCLUDE_PATH=$PREFIX/include \
 .

# Build step
ninja

# Install step
ninja install

# test
ninja test
