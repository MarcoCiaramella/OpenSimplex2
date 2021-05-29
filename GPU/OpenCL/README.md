# OpenSimplex2
GPU implementation in C of [OpenSimplex 2](https://github.com/KdotJPG/OpenSimplex2)

## How to use
Make sure you have a driver supporting OpenCL.
### OpenSimplex2F
```c
/* file test.c */

#include "OpenSimplex2F.h"

#define WIDTH 4096
#define HEIGHT 4096

int main(){
     OpenSimplexEnv* ose = initOpenSimplex();
     OpenSimplexGradients* osg = newOpenSimplexGradients(ose, 1234);
     OpenCLEnv openCLEnv = initOpenCL("OpenSimplex2F/OpenSimplex2F.cl", WIDTH, HEIGHT);

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
/* file test.c */

#include "OpenSimplex2S.h"

#define WIDTH 4096
#define HEIGHT 4096

int main(){
     OpenSimplexEnv* ose = initOpenSimplex();
     OpenSimplexGradients* osg = newOpenSimplexGradients(ose, 1234);
     OpenCLEnv openCLEnv = initOpenCL("OpenSimplex2S/OpenSimplex2S.cl", WIDTH, HEIGHT);

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
### gcc
```shell
gcc *.c /path/to/OpenSimplex2S/*.c -IC:/path/to/OpenCL-Headers -IC:/path/to/OpenSimplex2S -LC:/path/to/OpenCL_lib -lOpenCL -o test.exe
```