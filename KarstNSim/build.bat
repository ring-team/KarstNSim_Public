mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Debug -A x64 ..
cmake --build . --config Debug

cmake -DCMAKE_BUILD_TYPE=Release -A x64 ..
cmake --build . --config Release

pause