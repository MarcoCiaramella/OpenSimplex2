#include "..\opencl.h"



#define PSIZE 2048
#define PMASK 2047
#define N2 0.05481866495625118
#define N3 0.2781926117527186
#define N4 0.11127401889945551





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
    cl_short *perm;
    Grad2 *permGrad2;
    Grad3 *permGrad3;
    Grad4 *permGrad4;
} OpenSimplexGradients;

typedef struct {
    Grad2 *GRADIENTS_2D;
    Grad3 *GRADIENTS_3D;
    Grad4 *GRADIENTS_4D;
    LatticePoint2D* LOOKUP_2D;
    LatticePoint3D* LOOKUP_3D;
    LatticePoint4D* LOOKUP_4D;
} OpenSimplexEnv;

OpenSimplexEnv *initOpenSimplex();
OpenSimplexGradients *newOpenSimplexGradients(OpenSimplexEnv *ose, cl_long seed);
double *noise2(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise2_XBeforeY(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise3_Classic(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise3_XYBeforeZ(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise3_XZBeforeY(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise4_Classic(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise4_XYBeforeZW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise4_XZBeforeYW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
double *noise4_XYZBeforeW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg);
