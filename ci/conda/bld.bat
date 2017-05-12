mkdir build
cd build

REM Configure step
cmake -G "Ninja" -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
 -DCMAKE_BUILD_TYPE=RelWithDebInfo ^
 -DSMESH_EXTRA_WARNINGS=ON ^
 -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
 -DCMAKE_SYSTEM_PREFIX_PATH="%LIBRARY_PREFIX%" ^
 -DBoost_INCLUDE_DIRS="%LIBRARY_PREFIX%"\include/ ^
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
