# Leeds Loop Subdivision

## Build and Run
### QMake
To compile on feng-linux / feng-gps:
```bash
module add qt/5.13.0
qmake -project QT+=opengl
qmake
make
```
To run on feng-linux / feng-gps ( must has a GUI):
```bash
  ./LoopSubdivision  ./path_to/model.diredgenormal
```

### CMake
Cmake Support Qt6 upper version.
Input your own Qt Path in the `CMakeLists.txt` , in line `9`
```cmake
set(CMAKE_PREFIX_PATH "INPUT YOU OWN PATH")
```  

```bash
mkdir build
cd build
cmake ..
```  

### Test it
1. Press the "Loop Subdivision" Button to run the loop subdivision at one time. (on my machine, 1-4 times will get a correct answer, more than 5 times will take more time)
2. Press the "Save Current File" Button to write the current mesh into the hard disk. the write output file will exist in the same directory of the build directory.


