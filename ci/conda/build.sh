# Configure step
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=$PREFIX \
 -DCMAKE_BUILD_TYPE=RelWithDebInfo \
 -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG -D_DEBUG_ -Wall" \
 -DCMAKE_PREFIX_PATH=$PREFIX \
 -DCMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
 -DBoost_INCLUDE_DIRS=$PREFIX/include/ \
 -DSMESH_TESTING=ON \
 .

# Build step
ninja

# Install step
ninja install

# test
ninja test
