# OpenSimplex2
GPU implementation in C of [OpenSimplex 2](https://github.com/KdotJPG/OpenSimplex2)

## How to use
Make sure you have an OpenCL driver installed.
### OpenSimplex2F
```c
#include "OpenSimplex2F.h"

#define WIDTH 4096
#define HEIGHT 4096

int main(){
     OpenSimplexEnv* ose = initOpenSimplex();
     OpenSimplexGradients* osg = newOpenSimplexGradients(ose, 1234);
     OpenCLEnv openCLEnv = initOpenCL("OpenSimplex2F.cl", WIDTH, HEIGHT);

     double *output_buffer;

     output_buffer = noise2(&openCLEnv, ose, osg);
     output_buffer = noise2_XBeforeY(&openCLEnv, ose, osg);
     output_buffer = noise3_Classic(&openCLEnv, ose, osg);
     output_buffer = noise3_XYBeforeZ(&openCLEnv, ose, osg);
     output_buffer = noise3_XZBeforeY(&openCLEnv, ose, osg);
     output_buffer = noise4_Classic(&openCLEnv, ose, osg);
     output_buffer = noise4_XYBeforeZW(&openCLEnv, ose, osg);
     output_buffer = noise4_XZBeforeYW(&openCLEnv, ose, osg);
     output_buffer = noise4_XYZBeforeW(&openCLEnv, ose, osg);

     releaseOpenCL(&openCLEnv);

     return 0;
}
```
### OpenSimplex2S
```c
#include "OpenSimplex2S.h"

#define WIDTH 4096
#define HEIGHT 4096

int main(){
     OpenSimplexEnv* ose = initOpenSimplex();
     OpenSimplexGradients* osg = newOpenSimplexGradients(ose, 1234);
     OpenCLEnv openCLEnv = initOpenCL("OpenSimplex2S.cl", WIDTH, HEIGHT);

     double *output_buffer;

     output_buffer = noise2(&openCLEnv, ose, osg);
     output_buffer = noise2_XBeforeY(&openCLEnv, ose, osg);
     output_buffer = noise3_Classic(&openCLEnv, ose, osg);
     output_buffer = noise3_XYBeforeZ(&openCLEnv, ose, osg);
     output_buffer = noise3_XZBeforeY(&openCLEnv, ose, osg);
     output_buffer = noise4_Classic(&openCLEnv, ose, osg);
     output_buffer = noise4_XYBeforeZW(&openCLEnv, ose, osg);
     output_buffer = noise4_XZBeforeYW(&openCLEnv, ose, osg);
     output_buffer = noise4_XYZBeforeW(&openCLEnv, ose, osg);

     releaseOpenCL(&openCLEnv);
     
     return 0;
}
```

## How to compile
Download the official OpenCL headers [repository](https://github.com/KhronosGroup/OpenCL-Headers)
### gcc
Compiling OpenSimplex2F example in test directory:
```shell
cd GPU/OpenCL/test
gcc ../opencl.c ../OpenSimplex2F.c test_OpenSimplex2F.c bitmap.c -IC:/path/to/OpenCL-Headers -LC:/path/to/OpenCL_lib -lOpenCL -o test_OpenSimplex2F.exe
```
Compiling OpenSimplex2S example in test directory:
```shell
cd GPU/OpenCL/test
gcc ../opencl.c ../OpenSimplex2S.c test_OpenSimplex2S.c bitmap.c -IC:/path/to/OpenCL-Headers -LC:/path/to/OpenCL_lib -lOpenCL -o test_OpenSimplex2S.exe
```