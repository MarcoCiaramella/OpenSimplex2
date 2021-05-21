#include <CL/cl.h>


#define PSIZE 2048
#define PMASK 2047
#define N2 0.01001634121365712
#define N3 0.030485933181293584
#define N4 0.009202377986303158


#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD






typedef struct {
    cl_int xsv;
    cl_int ysv;
    cl_double dx;
    cl_double dy;
} LatticePoint2D;

typedef struct {
    cl_double dxr;
    cl_double dyr;
    cl_double dzr;
    cl_int xrv;
    cl_int yrv;
    cl_int zrv;
	cl_bool is_null;
} _LatticePoint3D;

typedef struct {
    _LatticePoint3D _this;
    _LatticePoint3D nextOnFailure;
    _LatticePoint3D nextOnSuccess;
} LatticePoint3D;

typedef struct {
    cl_int xsv;
    cl_int ysv;
    cl_int zsv;
    cl_int wsv;
    cl_double dx;
    cl_double dy;
    cl_double dz;
    cl_double dw;
    cl_double xsi;
    cl_double ysi;
    cl_double zsi;
    cl_double wsi;
    cl_double ssiDelta;
} LatticePoint4D;

typedef struct {
    cl_double dx;
    cl_double dy;
} Grad2;

typedef struct {
    cl_double dx;
    cl_double dy;
    cl_double dz;
} Grad3;

typedef struct {
    cl_double dx;
    cl_double dy;
    cl_double dz;
    cl_double dw;
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
    LatticePoint2D LOOKUP_2D[4];
    LatticePoint3D LOOKUP_3D[8];
    LatticePoint4D VERTICES_4D[16];
} OpenSimplexEnv;


OpenSimplexEnv initOpenSimplex();
OpenSimplexGradients newOpenSimplexGradients(OpenSimplexEnv *ose, cl_long seed);
