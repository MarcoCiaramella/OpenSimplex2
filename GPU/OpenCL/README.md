# OpenSimplex2
C implementation for GPU using OpenCL

## Note
All functions accept an input array of 2D, 3D and 4D points.
For example for noise2 you have to allocate and load an array of N points in the form of `[x1,y1,x2,y2,...,xn,yn]`.
With noise3 the input array must be in the form of `[x1,y1,z1,x2,y2,z2,...,xn,yn,zn]` and with noise4 `[x1,y1,z1,w1,x2,y2,z2,w2,...,xn,yn,zn,wn]`.
The ouput array stores the noise generated for each point.

## How to use
Make sure you have an OpenCL driver installed.
### OpenSimplex2F
You can find an [example](test/test_OpenSimplex2F.c) under test directory.
### OpenSimplex2S
You can find an [example](test/test_OpenSimplex2S.c) under test directory.

## How to compile
Download the official OpenCL headers [repository](https://github.com/KhronosGroup/OpenCL-Headers)
### gcc
Compiling OpenSimplex2F example:
```shell
cd GPU/OpenCL/test ; gcc ../opencl.c ../OpenSimplex2F.c test_OpenSimplex2F.c bitmap.c -IC:/path/to/OpenCL-Headers -LC:/path/to/OpenCL_lib -lOpenCL -o test_OpenSimplex2F.exe
```
Compiling OpenSimplex2S example:
```shell
cd GPU/OpenCL/test ; gcc ../opencl.c ../OpenSimplex2S.c test_OpenSimplex2S.c bitmap.c -IC:/path/to/OpenCL-Headers -LC:/path/to/OpenCL_lib -lOpenCL -o test_OpenSimplex2S.exe
```
