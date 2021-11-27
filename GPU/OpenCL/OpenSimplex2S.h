#include "opencl.h"



#define PSIZE 2048
#define N2 0.01001634121365712
#define N3 0.030485933181293584
#define N4 0.009202377986303158





typedef struct {
    cl_int xsv, ysv;
	cl_double dx, dy;
} LatticePoint2D;

typedef struct {
    cl_double dxr, dyr, dzr;
	cl_int xrv, yrv, zrv;
    cl_int nextOnFailure;
    cl_int nextOnSuccess;
} LatticePoint3D;

typedef struct {
    cl_int xsv, ysv, zsv, wsv;
	cl_double dx, dy, dz, dw;
} LatticePoint4D;

typedef struct {
    cl_double dx, dy;
} Grad2;

typedef struct {
    cl_double dx, dy, dz;
} Grad3;

typedef struct {
    cl_double dx, dy, dz, dw;
} Grad4;

typedef struct {
    cl_short perm[PSIZE];
    Grad2 permGrad2[PSIZE];
    Grad3 permGrad3[PSIZE];
    Grad4 permGrad4[PSIZE];
} OpenSimplexGradients;

typedef struct {
    Grad2 GRADIENTS_2D[PSIZE];
    Grad3 GRADIENTS_3D[PSIZE];
    Grad4 GRADIENTS_4D[PSIZE];
    LatticePoint2D LOOKUP_2D[32];
    LatticePoint3D LOOKUP_3D[112];
    LatticePoint4D LOOKUP_4D[3476];
} OpenSimplexEnv;

OpenSimplexEnv *initOpenSimplex();
OpenSimplexGradients *newOpenSimplexGradients(OpenSimplexEnv *ose, cl_long seed);
double *noise2(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise2_XBeforeY(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise3_Classic(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise3_XYBeforeZ(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise3_XZBeforeY(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise4_Classic(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise4_XYBeforeZW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise4_XZBeforeYW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
double *noise4_XYZBeforeW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg, double* input_buffer, size_t size_input_buffer);
