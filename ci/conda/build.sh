# Configure step
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=$PREFIX \
 -DCMAKE_BUILD_TYPE=RelWithDebInfo \
 -DRAISE_EXCEPTION_ON_FAILURE=ON \
 -DCMAKE_PREFIX_PATH=$PREFIX \
 -DCMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
 -DSMESH_USE_BUNDLED_BOOST=ON \
 -DSMESH_TESTING=ON \
 .

# Build step
ninja

# Install step
ninja install > installed_files.txt

# test
#ninja test
ctest --verbose
