#include <stdbool.h>


#define PSIZE 2048
#define PMASK 2047
#define N2 0.01001634121365712
#define N3 0.030485933181293584
#define N4 0.009202377986303158


#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD




#pragma pack()
typedef struct {
    int xsv, ysv;
    double dx, dy;
} LatticePoint2D;

typedef struct {
    double dxr, dyr, dzr;
    int xrv, yrv, zrv;
	bool is_null;
} _LatticePoint3D;

typedef struct {
    _LatticePoint3D _this;
    _LatticePoint3D nextOnFailure;
    _LatticePoint3D nextOnSuccess;
} LatticePoint3D;

typedef struct {
    int xsv, ysv, zsv, wsv;
    double dx, dy, dz, dw;
    double xsi, ysi, zsi, wsi;
    double ssiDelta;
} LatticePoint4D;

typedef struct {
    double dx, dy;
} Grad2;

typedef struct {
    double dx, dy, dz;
} Grad3;

typedef struct {
    double dx, dy, dz, dw;
} Grad4;

typedef struct {
    short perm[PSIZE];
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
OpenSimplexGradients newOpenSimplexGradients(OpenSimplexEnv *ose, long seed);
