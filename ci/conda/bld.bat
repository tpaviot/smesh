mkdir build
cd build

REM Configure step
cmake -G "Ninja" -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
 -DCMAKE_BUILD_TYPE=RelWithDebInfo ^
 -DRAISE_EXCEPTION_ON_FAILURE=ON ^
 -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
 -DCMAKE_SYSTEM_PREFIX_PATH="%LIBRARY_PREFIX%" ^
 -DSMESH_USE_BUNDLED_BOOST=ON ^
 -DOCE_INCLUDE_PATH="%LIBRARY_PREFIX%"\include\oce/ ^
 -DOCE_LIB_PATH="%LIBRARY_PREFIX%"\lib/ ^
 ..
if errorlevel 1 exit 1
 
REM Build step 
ninja
if errorlevel 1 exit 1

REM Install step
ninja install > installed_files.txt
if errorlevel 1 exit 1

REM Test step
ninja test
if errorlevel 1 exit 1
