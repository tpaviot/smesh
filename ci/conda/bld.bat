mkdir build
cd build

REM Configure step
cmake -G "Ninja" -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
 -DCMAKE_BUILD_TYPE=Release ^
 -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
 -DCMAKE_SYSTEM_PREFIX_PATH="%LIBRARY_PREFIX%" ^
 -DSMESH_BUILD_NETGENPLUGIN=ON ^
 -DNETGEN_NGLIB_INCLUDE_PATH="%LIBRARY_PREFIX%"\include\netgen ^
 -DBoost_INCLUDE_DIRS="%LIBRARY_PREFIX%"\include/ ^
 ..
if errorlevel 1 exit 1
 
REM Build step 
ninja
if errorlevel 1 exit 1

REM Install step
ninja install
if errorlevel 1 exit 1
