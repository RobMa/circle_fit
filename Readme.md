# circle_fit
**circle_fit** is a least-squares circle fitting library written in C++.

## Build instructions
Build requirements are a C++14 compiler, _cmake_ and _Python3_ for code generation.
```
git clone git@github.com:RobMa/circle_fit.git
cd circle_fit
mkdir build && cd build && cmake .. && cd ..
cmake --build build
```
The library depends on _Eigen3_, _spdlog_ and _Catch2_, _cmake_ will download these dependencies automatically.


