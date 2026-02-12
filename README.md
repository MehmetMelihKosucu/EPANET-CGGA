# EPANET-CGGA
Comprehensive Gradient Global Algorithm Extension of the EPANET

- If an adaptive hydraulic analysis is to be performed, the Hyd_Solver option must be specified as CGGA in the EPANET 3 Input file.

- The wave speed for pipes is determined using the getWaveSpeed() const function in the pipe.cpp file.


## Building
To build, use CMake on Windows with Visual Studio (tested with Visual Studio 2013 and validated with 2019):
```
mkdir build && cd build
cmake -G "Visual Studio n yyyy" ..
cmake --build . --config Release
```
## Contributors
```
Mehmet Melih Koşucu,		Istanbul Technical University

```

