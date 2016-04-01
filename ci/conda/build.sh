# Configure step
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_PREFIX_PATH=$PREFIX \
 -DCMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
 -DSMESH_TESTING=ON \
 .

# Build step
make -j 4

# test
if [ `uname` == Linux ]; then
  make test
fi


# Install step
# redirect stdout stream to a file because of many boost header files
make install > installed_files.txt

