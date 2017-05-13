# Configure step
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=$PREFIX \
 -DCMAKE_BUILD_TYPE=RelWithDebInfo \
 -DSMESH_EXTRA_WARNINGS=ON \
 -DCMAKE_PREFIX_PATH=$PREFIX \
 -DCMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
 -DBoost_INCLUDE_DIRS=$PREFIX/include/ \
 -DSMESH_TESTING=ON \
 .

# Build step
ninja

# Install step
ninja install > installed_files.txt

# test
ninja test
