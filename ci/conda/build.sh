# Configure step
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=$PREFIX \
 -DCMAKE_BUILD_TYPE=Release \
 -DRAISE_EXCEPTION_ON_FAILURE=ON \
 -DCMAKE_PREFIX_PATH=$PREFIX \
 -DCMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
 -DSMESH_TESTING=ON \
 .

# Build step
ninja

# Install step
ninja install > installed_files.txt

# test
ninja test
