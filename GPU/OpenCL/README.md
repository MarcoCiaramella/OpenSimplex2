# OpenSimplex2
C implementation for GPU using OpenCL

## Note
The library in this early version generates a 2D noise over a grid WIDTH x HEIGHT so you can use the result as image only (for example as a heightmap).

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
     OpenCLEnv openCLEnv = loadOpenCL("OpenSimplex2F.cl", WIDTH, HEIGHT);

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
     OpenCLEnv openCLEnv = loadOpenCL("OpenSimplex2S.cl", WIDTH, HEIGHT);

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
Compiling OpenSimplex2F example:
```shell
cd GPU/OpenCL/test ; gcc ../opencl.c ../OpenSimplex2F.c test_OpenSimplex2F.c bitmap.c -IC:/path/to/OpenCL-Headers -LC:/path/to/OpenCL_lib -lOpenCL -o test_OpenSimplex2F.exe
```
Compiling OpenSimplex2S example:
```shell
cd GPU/OpenCL/test ; gcc ../opencl.c ../OpenSimplex2S.c test_OpenSimplex2S.c bitmap.c -IC:/path/to/OpenCL-Headers -LC:/path/to/OpenCL_lib -lOpenCL -o test_OpenSimplex2S.exe
```

## Performance
Image 4096x4096 on GPU NVIDIA GeForce GTX 1650
### OpenSimplex2F
* noise2 0.070000s
* noise2_XBeforeY 0.073000s
* noise3_Classic 0.094000s
* noise3_XYBeforeZ 0.095000s
* noise3_XZBeforeY 0.101000s
* noise4_Classic 0.154000s
* noise4_XYBeforeZW 0.125000s
* noise4_XZBeforeYW 0.120000s
* noise4_XYZBeforeW 0.121000s
### OpenSimplex2S
* noise2 0.082000s
* noise2_XBeforeY 0.080000s
* noise3_Classic 0.112000s
* noise3_XYBeforeZ 0.110000s
* noise3_XZBeforeY 0.115000s
* noise4_Classic 0.175000s
* noise4_XYBeforeZW 0.143000s
* noise4_XZBeforeYW 0.136000s
* noise4_XYZBeforeW 0.150000s
