#include "OpenSimplex2S.h"
#include <stdbool.h>






Grad2 *_newGrad2Arr(unsigned int size){
    return (Grad2 *) malloc(sizeof(Grad2)*size);
}

Grad3 *_newGrad3Arr(unsigned int size){
    return (Grad3 *) malloc(sizeof(Grad3)*size);
}

Grad4 *_newGrad4Arr(unsigned int size){
    return (Grad4 *) malloc(sizeof(Grad4)*size);
}

short *_newShortArr(unsigned int size){
    return (short *) malloc(sizeof(short)*size);
}

LatticePoint2D _newLatticePoint2D(int xsv, int ysv){
	LatticePoint2D lp2D;
	lp2D.xsv = xsv;
    lp2D.ysv = ysv;
	double ssv = (xsv + ysv) * -0.211324865405187;
	lp2D.dx = -xsv - ssv;
	lp2D.dy = -ysv - ssv;
	return lp2D;
}

LatticePoint3D _newLatticePoint3D(int xrv, int yrv, int zrv, int lattice){
	LatticePoint3D lp3D;
	lp3D.dxr = -xrv + lattice * 0.5;
    lp3D.dyr = -yrv + lattice * 0.5;
    lp3D.dzr = -zrv + lattice * 0.5;
	lp3D.xrv = xrv + lattice * 1024;
    lp3D.yrv = yrv + lattice * 1024;
    lp3D.zrv = zrv + lattice * 1024;
	return lp3D;
}

LatticePoint4D _newLatticePoint4D(int xsv, int ysv, int zsv, int wsv){
	LatticePoint4D lp4D;
    lp4D.xsv = xsv;
    lp4D.ysv = ysv;
    lp4D.zsv = zsv;
    lp4D.wsv = wsv;
	double ssv = (xsv + ysv + zsv + wsv) * -0.138196601125011;
	lp4D.dx = -xsv - ssv;
	lp4D.dy = -ysv - ssv;
	lp4D.dz = -zsv - ssv;
	lp4D.dw = -wsv - ssv;
	return lp4D;
}

LatticePoint2D* _newLatticePoint2DConstArray(){
	LatticePoint2D* arr = (LatticePoint2D*) malloc(sizeof(LatticePoint2D) * 8 * 4);
	for (int i = 0; i < 8; i++){
		int i1, j1, i2, j2;
		if ((i & 1) == 0){
			if ((i & 2) == 0){
				i1 = -1;
				j1 = 0;
			}
			else{
				i1 = 1;
				j1 = 0;
			}
			if ((i & 4) == 0){
				i2 = 0;
				j2 = -1;
			}
			else{
				i2 = 0;
				j2 = 1;
			}
		}
		else{
			if ((i & 2) != 0){
				i1 = 2;
				j1 = 1;
			}
			else{
				i1 = 0;
				j1 = 1;
			}
			if ((i & 4) != 0){
				i2 = 1;
				j2 = 2;
			}
			else{
				i2 = 1;
				j2 = 0;
			}
		}
		arr[i * 4 + 0] = _newLatticePoint2D(0, 0);
		arr[i * 4 + 1] = _newLatticePoint2D(1, 1);
		arr[i * 4 + 2] = _newLatticePoint2D(i1, j1);
		arr[i * 4 + 3] = _newLatticePoint2D(i2, j2);
	}
	return arr;
}

LatticePoint3D* _newLatticePoint3DConstArray(){
	LatticePoint3D* arr = (LatticePoint3D*) malloc(sizeof(LatticePoint3D) * 8 * 13);
	int j = 7;
	for (int i = 0; i < 8; i++){
		int i1, j1, k1, i2, j2, k2;
		i1 = (i >> 0) & 1;
		j1 = (i >> 1) & 1;
		k1 = (i >> 2) & 1;
		i2 = i1 ^ 1;
		j2 = j1 ^ 1;
		k2 = k1 ^ 1;

		// The two points within this octant, one from each of the two cubic half-lattices.
		LatticePoint3D c0 = _newLatticePoint3D(i1, j1, k1, 0);
		LatticePoint3D c1 = _newLatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1);

		// (1, 0, 0) vs (0, 1, 1) away from octant.
		LatticePoint3D c2 = _newLatticePoint3D(i1 ^ 1, j1, k1, 0);
		LatticePoint3D c3 = _newLatticePoint3D(i1, j1 ^ 1, k1 ^ 1, 0);

		// (1, 0, 0) vs (0, 1, 1) away from octant, on second half-lattice.
		LatticePoint3D c4 = _newLatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1);
		LatticePoint3D c5 = _newLatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + (k2 ^ 1), 1);

		// (0, 1, 0) vs (1, 0, 1) away from octant.
		LatticePoint3D c6 = _newLatticePoint3D(i1, j1 ^ 1, k1, 0);
		LatticePoint3D c7 = _newLatticePoint3D(i1 ^ 1, j1, k1 ^ 1, 0);

		// (0, 1, 0) vs (1, 0, 1) away from octant, on second half-lattice.
		LatticePoint3D c8 = _newLatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1);
		LatticePoint3D c9 = _newLatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + (k2 ^ 1), 1);

		// (0, 0, 1) vs (1, 1, 0) away from octant.
		LatticePoint3D cA = _newLatticePoint3D(i1, j1, k1 ^ 1, 0);
		LatticePoint3D cB = _newLatticePoint3D(i1 ^ 1, j1 ^ 1, k1, 0);

		// (0, 0, 1) vs (1, 1, 0) away from octant, on second half-lattice.
		LatticePoint3D cC = _newLatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1);
		LatticePoint3D cD = _newLatticePoint3D(i1 + (i2 ^ 1), j1 + (j2 ^ 1), k1 + k2, 1);

		// First two points are guaranteed.
		c0.nextOnFailure = c0.nextOnSuccess = j+1;
		c1.nextOnFailure = c1.nextOnSuccess = j+2;

		// If c2 is in range, then we know c3 and c4 are not.
		c2.nextOnFailure = j+3;
		c2.nextOnSuccess = j+5;
		c3.nextOnFailure = j+4;
		c3.nextOnSuccess = j+4;

		// If c4 is in range, then we know c5 is not.
		c4.nextOnFailure = j+5;
		c4.nextOnSuccess = j+6;
		c5.nextOnFailure = c5.nextOnSuccess = j+6;

		// If c6 is in range, then we know c7 and c8 are not.
		c6.nextOnFailure = j+7;
		c6.nextOnSuccess = j+9;
		c7.nextOnFailure = j+8;
		c7.nextOnSuccess = j+8;

		// If c8 is in range, then we know c9 is not.
		c8.nextOnFailure = j+9;
		c8.nextOnSuccess = j+10;
		c9.nextOnFailure = c9.nextOnSuccess = j+10;

		// If cA is in range, then we know cB and cC are not.
		cA.nextOnFailure = j+11;
		cA.nextOnSuccess = j+13;
		cB.nextOnFailure = j+12;
		cB.nextOnSuccess = j+12;

		// If cC is in range, then we know cD is not.
		cC.nextOnFailure = j+13;
		cC.nextOnSuccess = -1;
		cD.nextOnFailure = cD.nextOnSuccess = -1;

		arr[i] = c0;
		arr[++j] = c1;
		arr[++j] = c2;
		arr[++j] = c3;
		arr[++j] = c4;
		arr[++j] = c5;
		arr[++j] = c6;
		arr[++j] = c7;
		arr[++j] = c8;
		arr[++j] = c9;
		arr[++j] = cA;
		arr[++j] = cB;
		arr[++j] = cC;
		arr[++j] = cD;
	}
	return arr;
}

vect _newVect(unsigned int length, int *data){
	vect v = {length, data};
	return v;
}

vect *_newLatticePoint4DConstArray(){
	vect *plp4DArr = (vect *) malloc(sizeof(vect) * 256);
	vect lookup4DPregen[256];
	int i = 0;
	{
		int arr[] = {0x15, 0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x01, 0x05, 0x11, 0x15, 0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x01, 0x15, 0x16, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x04, 0x05, 0x14, 0x15, 0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x6A, 0x9A, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x04, 0x15, 0x19, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x5E, 0x6A, 0x9A, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x05, 0x15, 0x1A, 0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x5E, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x16, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x6B, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x6B, 0x9A, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x19, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x9A, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x1A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x10, 0x11, 0x14, 0x15, 0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x67, 0x6A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x16, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x6B, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x6D, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x19, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x6E, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x10, 0x15, 0x25, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0x76, 0xA6, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x11, 0x15, 0x26, 0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x67, 0x6A, 0x76, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x25, 0x55, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x25, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA6, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x26, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0x79, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x25, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x25, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x14, 0x15, 0x29, 0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0x6D, 0x79, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x29, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x7A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x46, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x9A, 0x9B, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x49, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x99, 0x9A, 0x9E, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x59, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x56, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xA7, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAD, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x6A, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x46, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x95, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x49, 0x55, 0x59, 0x5A, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x52, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x58, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xAF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x61, 0x65, 0x66, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x64, 0x65, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xBB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x40, 0x45, 0x51, 0x54, 0x55, 0x85, 0x91, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xD6, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x86, 0x92, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB, 0xD6, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xDA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xDA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x86, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xD9, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x59, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xDA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xDA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x89, 0x95, 0x98, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE, 0xD9, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x59, 0x89, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xDA, 0xEA, 0xEF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x91, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x92, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x94, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x95, 0x98, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xEA, 0xEF};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xE5, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x65, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x65, 0x94, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x94, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA, 0xEB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA1, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA, 0xE5, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x65, 0x95, 0xA1, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xE6, 0xEA, 0xFB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x65, 0x95, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xEA, 0xFB};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xE9, 0xEA, 0xFE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA, 0xFE};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	{
		int arr[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA, 0xEA};
		lookup4DPregen[i++] = _newVect(sizeof arr / sizeof *arr, arr);
	}
	LatticePoint4D* latticePoints = (LatticePoint4D*) malloc(sizeof(LatticePoint4D) * 256);
	for (int i = 0; i < 256; i++){
		int cx = ((i >> 0) & 3) - 1;
		int cy = ((i >> 2) & 3) - 1;
		int cz = ((i >> 4) & 3) - 1;
		int cw = ((i >> 6) & 3) - 1;
		latticePoints[i] = _newLatticePoint4D(cx, cy, cz, cw);
	}
	for (int i = 0; i < 256; i++){
		vect v;
		v.data = (LatticePoint4D *)malloc(sizeof(LatticePoint4D) * lookup4DPregen[i].length);
		v.length = lookup4DPregen[i].length;
		for (int j = 0; j < lookup4DPregen[i].length; j++){
			((LatticePoint4D *)v.data)[j] = *(latticePoints[((int *)lookup4DPregen[i].data)[j]]);
		}
		plp4DArr[i] = v;
	}
	return plp4DArr;
}

Grad2 _newGrad2(double dx, double dy){
    Grad2 grad2;
    grad2.dx = dx;
    grad2.dy = dy;
    return grad2;
}

Grad3 _newGrad3(double dx, double dy, double dz){
    Grad3 grad3;
    grad3.dx = dx;
    grad3.dy = dy;
    grad3.dz = dz;
    return grad3;
}

Grad4 _newGrad4(double dx, double dy, double dz, double dw){
    Grad4 grad4;
    grad4.dx = dx;
    grad4.dy = dy;
    grad4.dz = dz;
    grad4.dw = dw;
    return grad4;
}

Grad2 *_newGrad2ConstArray(){
    Grad2 *arr = (Grad2 *) malloc(sizeof(Grad2)*24);
    int i = 0;
	arr[i++] = _newGrad2(0.130526192220052, 0.99144486137381);
	arr[i++] = _newGrad2(0.38268343236509, 0.923879532511287);
	arr[i++] = _newGrad2(0.608761429008721, 0.793353340291235);
	arr[i++] = _newGrad2(0.793353340291235, 0.608761429008721);
	arr[i++] = _newGrad2(0.923879532511287, 0.38268343236509);
	arr[i++] = _newGrad2(0.99144486137381, 0.130526192220051);
	arr[i++] = _newGrad2(0.99144486137381, -0.130526192220051);
	arr[i++] = _newGrad2(0.923879532511287, -0.38268343236509);
	arr[i++] = _newGrad2(0.793353340291235, -0.60876142900872);
	arr[i++] = _newGrad2(0.608761429008721, -0.793353340291235);
	arr[i++] = _newGrad2(0.38268343236509, -0.923879532511287);
	arr[i++] = _newGrad2(0.130526192220052, -0.99144486137381);
	arr[i++] = _newGrad2(-0.130526192220052, -0.99144486137381);
	arr[i++] = _newGrad2(-0.38268343236509, -0.923879532511287);
	arr[i++] = _newGrad2(-0.608761429008721, -0.793353340291235);
	arr[i++] = _newGrad2(-0.793353340291235, -0.608761429008721);
	arr[i++] = _newGrad2(-0.923879532511287, -0.38268343236509);
	arr[i++] = _newGrad2(-0.99144486137381, -0.130526192220052);
	arr[i++] = _newGrad2(-0.99144486137381, 0.130526192220051);
	arr[i++] = _newGrad2(-0.923879532511287, 0.38268343236509);
	arr[i++] = _newGrad2(-0.793353340291235, 0.608761429008721);
	arr[i++] = _newGrad2(-0.608761429008721, 0.793353340291235);
	arr[i++] = _newGrad2(-0.38268343236509, 0.923879532511287);
	arr[i++] = _newGrad2(-0.130526192220052, 0.99144486137381);
	Grad2 *gradients2D = _newGrad2Arr(PSIZE);
	for (int i = 0; i < 24; i++){
		arr[i].dx /= N2;
		arr[i].dy /= N2;
	}
	for (int i = 0; i < PSIZE; i++){
		gradients2D[i] = arr[i % 24];
	}
	return gradients2D;
}

Grad3 *_newGrad3ConstArray(){
	Grad3 *arr = (Grad3 *)malloc(sizeof(Grad3) * 48);
	int i = 0;
	arr[i++] = _newGrad3(-2.22474487139, -2.22474487139, -1.0);
	arr[i++] = _newGrad3(-2.22474487139, -2.22474487139, 1.0);
	arr[i++] = _newGrad3(-3.0862664687972017, -1.1721513422464978, 0.0);
	arr[i++] = _newGrad3(-1.1721513422464978, -3.0862664687972017, 0.0);
	arr[i++] = _newGrad3(-2.22474487139, -1.0, -2.22474487139);
	arr[i++] = _newGrad3(-2.22474487139, 1.0, -2.22474487139);
	arr[i++] = _newGrad3(-1.1721513422464978, 0.0, -3.0862664687972017);
	arr[i++] = _newGrad3(-3.0862664687972017, 0.0, -1.1721513422464978);
	arr[i++] = _newGrad3(-2.22474487139, -1.0, 2.22474487139);
	arr[i++] = _newGrad3(-2.22474487139, 1.0, 2.22474487139);
	arr[i++] = _newGrad3(-3.0862664687972017, 0.0, 1.1721513422464978);
	arr[i++] = _newGrad3(-1.1721513422464978, 0.0, 3.0862664687972017);
	arr[i++] = _newGrad3(-2.22474487139, 2.22474487139, -1.0);
	arr[i++] = _newGrad3(-2.22474487139, 2.22474487139, 1.0);
	arr[i++] = _newGrad3(-1.1721513422464978, 3.0862664687972017, 0.0);
	arr[i++] = _newGrad3(-3.0862664687972017, 1.1721513422464978, 0.0);
	arr[i++] = _newGrad3(-1.0, -2.22474487139, -2.22474487139);
	arr[i++] = _newGrad3(1.0, -2.22474487139, -2.22474487139);
	arr[i++] = _newGrad3(0.0, -3.0862664687972017, -1.1721513422464978);
	arr[i++] = _newGrad3(0.0, -1.1721513422464978, -3.0862664687972017);
	arr[i++] = _newGrad3(-1.0, -2.22474487139, 2.22474487139);
	arr[i++] = _newGrad3(1.0, -2.22474487139, 2.22474487139);
	arr[i++] = _newGrad3(0.0, -1.1721513422464978, 3.0862664687972017);
	arr[i++] = _newGrad3(0.0, -3.0862664687972017, 1.1721513422464978);
	arr[i++] = _newGrad3(-1.0, 2.22474487139, -2.22474487139);
	arr[i++] = _newGrad3(1.0, 2.22474487139, -2.22474487139);
	arr[i++] = _newGrad3(0.0, 1.1721513422464978, -3.0862664687972017);
	arr[i++] = _newGrad3(0.0, 3.0862664687972017, -1.1721513422464978);
	arr[i++] = _newGrad3(-1.0, 2.22474487139, 2.22474487139);
	arr[i++] = _newGrad3(1.0, 2.22474487139, 2.22474487139);
	arr[i++] = _newGrad3(0.0, 3.0862664687972017, 1.1721513422464978);
	arr[i++] = _newGrad3(0.0, 1.1721513422464978, 3.0862664687972017);
	arr[i++] = _newGrad3(2.22474487139, -2.22474487139, -1.0);
	arr[i++] = _newGrad3(2.22474487139, -2.22474487139, 1.0);
	arr[i++] = _newGrad3(1.1721513422464978, -3.0862664687972017, 0.0);
	arr[i++] = _newGrad3(3.0862664687972017, -1.1721513422464978, 0.0);
	arr[i++] = _newGrad3(2.22474487139, -1.0, -2.22474487139);
	arr[i++] = _newGrad3(2.22474487139, 1.0, -2.22474487139);
	arr[i++] = _newGrad3(3.0862664687972017, 0.0, -1.1721513422464978);
	arr[i++] = _newGrad3(1.1721513422464978, 0.0, -3.0862664687972017);
	arr[i++] = _newGrad3(2.22474487139, -1.0, 2.22474487139);
	arr[i++] = _newGrad3(2.22474487139, 1.0, 2.22474487139);
	arr[i++] = _newGrad3(1.1721513422464978, 0.0, 3.0862664687972017);
	arr[i++] = _newGrad3(3.0862664687972017, 0.0, 1.1721513422464978);
	arr[i++] = _newGrad3(2.22474487139, 2.22474487139, -1.0);
	arr[i++] = _newGrad3(2.22474487139, 2.22474487139, 1.0);
	arr[i++] = _newGrad3(3.0862664687972017, 1.1721513422464978, 0.0);
	arr[i++] = _newGrad3(1.1721513422464978, 3.0862664687972017, 0.0);
	Grad3 *gradients3D = _newGrad3Arr(PSIZE);
	for (int i = 0; i < 48; i++){
		arr[i].dx /= N3;
		arr[i].dy /= N3;
		arr[i].dz /= N3;
	}
	for (int i = 0; i < PSIZE; i++){
		gradients3D[i] = arr[i % 48];
	}
	return gradients3D;
}

Grad4 *_newGrad4ConstArray(){
	Grad4 *arr = (Grad4 *)malloc(sizeof(Grad4) * 160);
	int i = 0;
	arr[i++] = _newGrad4(-0.753341017856078, -0.37968289875261624, -0.37968289875261624, -0.37968289875261624);
	arr[i++] = _newGrad4(-0.7821684431180708, -0.4321472685365301, -0.4321472685365301, 0.12128480194602098);
	arr[i++] = _newGrad4(-0.7821684431180708, -0.4321472685365301, 0.12128480194602098, -0.4321472685365301);
	arr[i++] = _newGrad4(-0.7821684431180708, 0.12128480194602098, -0.4321472685365301, -0.4321472685365301);
	arr[i++] = _newGrad4(-0.8586508742123365, -0.508629699630796, 0.044802370851755174, 0.044802370851755174);
	arr[i++] = _newGrad4(-0.8586508742123365, 0.044802370851755174, -0.508629699630796, 0.044802370851755174);
	arr[i++] = _newGrad4(-0.8586508742123365, 0.044802370851755174, 0.044802370851755174, -0.508629699630796);
	arr[i++] = _newGrad4(-0.9982828964265062, -0.03381941603233842, -0.03381941603233842, -0.03381941603233842);
	arr[i++] = _newGrad4(-0.37968289875261624, -0.753341017856078, -0.37968289875261624, -0.37968289875261624);
	arr[i++] = _newGrad4(-0.4321472685365301, -0.7821684431180708, -0.4321472685365301, 0.12128480194602098);
	arr[i++] = _newGrad4(-0.4321472685365301, -0.7821684431180708, 0.12128480194602098, -0.4321472685365301);
	arr[i++] = _newGrad4(0.12128480194602098, -0.7821684431180708, -0.4321472685365301, -0.4321472685365301);
	arr[i++] = _newGrad4(-0.508629699630796, -0.8586508742123365, 0.044802370851755174, 0.044802370851755174);
	arr[i++] = _newGrad4(0.044802370851755174, -0.8586508742123365, -0.508629699630796, 0.044802370851755174);
	arr[i++] = _newGrad4(0.044802370851755174, -0.8586508742123365, 0.044802370851755174, -0.508629699630796);
	arr[i++] = _newGrad4(-0.03381941603233842, -0.9982828964265062, -0.03381941603233842, -0.03381941603233842);
	arr[i++] = _newGrad4(-0.37968289875261624, -0.37968289875261624, -0.753341017856078, -0.37968289875261624);
	arr[i++] = _newGrad4(-0.4321472685365301, -0.4321472685365301, -0.7821684431180708, 0.12128480194602098);
	arr[i++] = _newGrad4(-0.4321472685365301, 0.12128480194602098, -0.7821684431180708, -0.4321472685365301);
	arr[i++] = _newGrad4(0.12128480194602098, -0.4321472685365301, -0.7821684431180708, -0.4321472685365301);
	arr[i++] = _newGrad4(-0.508629699630796, 0.044802370851755174, -0.8586508742123365, 0.044802370851755174);
	arr[i++] = _newGrad4(0.044802370851755174, -0.508629699630796, -0.8586508742123365, 0.044802370851755174);
	arr[i++] = _newGrad4(0.044802370851755174, 0.044802370851755174, -0.8586508742123365, -0.508629699630796);
	arr[i++] = _newGrad4(-0.03381941603233842, -0.03381941603233842, -0.9982828964265062, -0.03381941603233842);
	arr[i++] = _newGrad4(-0.37968289875261624, -0.37968289875261624, -0.37968289875261624, -0.753341017856078);
	arr[i++] = _newGrad4(-0.4321472685365301, -0.4321472685365301, 0.12128480194602098, -0.7821684431180708);
	arr[i++] = _newGrad4(-0.4321472685365301, 0.12128480194602098, -0.4321472685365301, -0.7821684431180708);
	arr[i++] = _newGrad4(0.12128480194602098, -0.4321472685365301, -0.4321472685365301, -0.7821684431180708);
	arr[i++] = _newGrad4(-0.508629699630796, 0.044802370851755174, 0.044802370851755174, -0.8586508742123365);
	arr[i++] = _newGrad4(0.044802370851755174, -0.508629699630796, 0.044802370851755174, -0.8586508742123365);
	arr[i++] = _newGrad4(0.044802370851755174, 0.044802370851755174, -0.508629699630796, -0.8586508742123365);
	arr[i++] = _newGrad4(-0.03381941603233842, -0.03381941603233842, -0.03381941603233842, -0.9982828964265062);
	arr[i++] = _newGrad4(-0.6740059517812944, -0.3239847771997537, -0.3239847771997537, 0.5794684678643381);
	arr[i++] = _newGrad4(-0.7504883828755602, -0.4004672082940195, 0.15296486218853164, 0.5029860367700724);
	arr[i++] = _newGrad4(-0.7504883828755602, 0.15296486218853164, -0.4004672082940195, 0.5029860367700724);
	arr[i++] = _newGrad4(-0.8828161875373585, 0.08164729285680945, 0.08164729285680945, 0.4553054119602712);
	arr[i++] = _newGrad4(-0.4553054119602712, -0.08164729285680945, -0.08164729285680945, 0.8828161875373585);
	arr[i++] = _newGrad4(-0.5029860367700724, -0.15296486218853164, 0.4004672082940195, 0.7504883828755602);
	arr[i++] = _newGrad4(-0.5029860367700724, 0.4004672082940195, -0.15296486218853164, 0.7504883828755602);
	arr[i++] = _newGrad4(-0.5794684678643381, 0.3239847771997537, 0.3239847771997537, 0.6740059517812944);
	arr[i++] = _newGrad4(-0.3239847771997537, -0.6740059517812944, -0.3239847771997537, 0.5794684678643381);
	arr[i++] = _newGrad4(-0.4004672082940195, -0.7504883828755602, 0.15296486218853164, 0.5029860367700724);
	arr[i++] = _newGrad4(0.15296486218853164, -0.7504883828755602, -0.4004672082940195, 0.5029860367700724);
	arr[i++] = _newGrad4(0.08164729285680945, -0.8828161875373585, 0.08164729285680945, 0.4553054119602712);
	arr[i++] = _newGrad4(-0.08164729285680945, -0.4553054119602712, -0.08164729285680945, 0.8828161875373585);
	arr[i++] = _newGrad4(-0.15296486218853164, -0.5029860367700724, 0.4004672082940195, 0.7504883828755602);
	arr[i++] = _newGrad4(0.4004672082940195, -0.5029860367700724, -0.15296486218853164, 0.7504883828755602);
	arr[i++] = _newGrad4(0.3239847771997537, -0.5794684678643381, 0.3239847771997537, 0.6740059517812944);
	arr[i++] = _newGrad4(-0.3239847771997537, -0.3239847771997537, -0.6740059517812944, 0.5794684678643381);
	arr[i++] = _newGrad4(-0.4004672082940195, 0.15296486218853164, -0.7504883828755602, 0.5029860367700724);
	arr[i++] = _newGrad4(0.15296486218853164, -0.4004672082940195, -0.7504883828755602, 0.5029860367700724);
	arr[i++] = _newGrad4(0.08164729285680945, 0.08164729285680945, -0.8828161875373585, 0.4553054119602712);
	arr[i++] = _newGrad4(-0.08164729285680945, -0.08164729285680945, -0.4553054119602712, 0.8828161875373585);
	arr[i++] = _newGrad4(-0.15296486218853164, 0.4004672082940195, -0.5029860367700724, 0.7504883828755602);
	arr[i++] = _newGrad4(0.4004672082940195, -0.15296486218853164, -0.5029860367700724, 0.7504883828755602);
	arr[i++] = _newGrad4(0.3239847771997537, 0.3239847771997537, -0.5794684678643381, 0.6740059517812944);
	arr[i++] = _newGrad4(-0.6740059517812944, -0.3239847771997537, 0.5794684678643381, -0.3239847771997537);
	arr[i++] = _newGrad4(-0.7504883828755602, -0.4004672082940195, 0.5029860367700724, 0.15296486218853164);
	arr[i++] = _newGrad4(-0.7504883828755602, 0.15296486218853164, 0.5029860367700724, -0.4004672082940195);
	arr[i++] = _newGrad4(-0.8828161875373585, 0.08164729285680945, 0.4553054119602712, 0.08164729285680945);
	arr[i++] = _newGrad4(-0.4553054119602712, -0.08164729285680945, 0.8828161875373585, -0.08164729285680945);
	arr[i++] = _newGrad4(-0.5029860367700724, -0.15296486218853164, 0.7504883828755602, 0.4004672082940195);
	arr[i++] = _newGrad4(-0.5029860367700724, 0.4004672082940195, 0.7504883828755602, -0.15296486218853164);
	arr[i++] = _newGrad4(-0.5794684678643381, 0.3239847771997537, 0.6740059517812944, 0.3239847771997537);
	arr[i++] = _newGrad4(-0.3239847771997537, -0.6740059517812944, 0.5794684678643381, -0.3239847771997537);
	arr[i++] = _newGrad4(-0.4004672082940195, -0.7504883828755602, 0.5029860367700724, 0.15296486218853164);
	arr[i++] = _newGrad4(0.15296486218853164, -0.7504883828755602, 0.5029860367700724, -0.4004672082940195);
	arr[i++] = _newGrad4(0.08164729285680945, -0.8828161875373585, 0.4553054119602712, 0.08164729285680945);
	arr[i++] = _newGrad4(-0.08164729285680945, -0.4553054119602712, 0.8828161875373585, -0.08164729285680945);
	arr[i++] = _newGrad4(-0.15296486218853164, -0.5029860367700724, 0.7504883828755602, 0.4004672082940195);
	arr[i++] = _newGrad4(0.4004672082940195, -0.5029860367700724, 0.7504883828755602, -0.15296486218853164);
	arr[i++] = _newGrad4(0.3239847771997537, -0.5794684678643381, 0.6740059517812944, 0.3239847771997537);
	arr[i++] = _newGrad4(-0.3239847771997537, -0.3239847771997537, 0.5794684678643381, -0.6740059517812944);
	arr[i++] = _newGrad4(-0.4004672082940195, 0.15296486218853164, 0.5029860367700724, -0.7504883828755602);
	arr[i++] = _newGrad4(0.15296486218853164, -0.4004672082940195, 0.5029860367700724, -0.7504883828755602);
	arr[i++] = _newGrad4(0.08164729285680945, 0.08164729285680945, 0.4553054119602712, -0.8828161875373585);
	arr[i++] = _newGrad4(-0.08164729285680945, -0.08164729285680945, 0.8828161875373585, -0.4553054119602712);
	arr[i++] = _newGrad4(-0.15296486218853164, 0.4004672082940195, 0.7504883828755602, -0.5029860367700724);
	arr[i++] = _newGrad4(0.4004672082940195, -0.15296486218853164, 0.7504883828755602, -0.5029860367700724);
	arr[i++] = _newGrad4(0.3239847771997537, 0.3239847771997537, 0.6740059517812944, -0.5794684678643381);
	arr[i++] = _newGrad4(-0.6740059517812944, 0.5794684678643381, -0.3239847771997537, -0.3239847771997537);
	arr[i++] = _newGrad4(-0.7504883828755602, 0.5029860367700724, -0.4004672082940195, 0.15296486218853164);
	arr[i++] = _newGrad4(-0.7504883828755602, 0.5029860367700724, 0.15296486218853164, -0.4004672082940195);
	arr[i++] = _newGrad4(-0.8828161875373585, 0.4553054119602712, 0.08164729285680945, 0.08164729285680945);
	arr[i++] = _newGrad4(-0.4553054119602712, 0.8828161875373585, -0.08164729285680945, -0.08164729285680945);
	arr[i++] = _newGrad4(-0.5029860367700724, 0.7504883828755602, -0.15296486218853164, 0.4004672082940195);
	arr[i++] = _newGrad4(-0.5029860367700724, 0.7504883828755602, 0.4004672082940195, -0.15296486218853164);
	arr[i++] = _newGrad4(-0.5794684678643381, 0.6740059517812944, 0.3239847771997537, 0.3239847771997537);
	arr[i++] = _newGrad4(-0.3239847771997537, 0.5794684678643381, -0.6740059517812944, -0.3239847771997537);
	arr[i++] = _newGrad4(-0.4004672082940195, 0.5029860367700724, -0.7504883828755602, 0.15296486218853164);
	arr[i++] = _newGrad4(0.15296486218853164, 0.5029860367700724, -0.7504883828755602, -0.4004672082940195);
	arr[i++] = _newGrad4(0.08164729285680945, 0.4553054119602712, -0.8828161875373585, 0.08164729285680945);
	arr[i++] = _newGrad4(-0.08164729285680945, 0.8828161875373585, -0.4553054119602712, -0.08164729285680945);
	arr[i++] = _newGrad4(-0.15296486218853164, 0.7504883828755602, -0.5029860367700724, 0.4004672082940195);
	arr[i++] = _newGrad4(0.4004672082940195, 0.7504883828755602, -0.5029860367700724, -0.15296486218853164);
	arr[i++] = _newGrad4(0.3239847771997537, 0.6740059517812944, -0.5794684678643381, 0.3239847771997537);
	arr[i++] = _newGrad4(-0.3239847771997537, 0.5794684678643381, -0.3239847771997537, -0.6740059517812944);
	arr[i++] = _newGrad4(-0.4004672082940195, 0.5029860367700724, 0.15296486218853164, -0.7504883828755602);
	arr[i++] = _newGrad4(0.15296486218853164, 0.5029860367700724, -0.4004672082940195, -0.7504883828755602);
	arr[i++] = _newGrad4(0.08164729285680945, 0.4553054119602712, 0.08164729285680945, -0.8828161875373585);
	arr[i++] = _newGrad4(-0.08164729285680945, 0.8828161875373585, -0.08164729285680945, -0.4553054119602712);
	arr[i++] = _newGrad4(-0.15296486218853164, 0.7504883828755602, 0.4004672082940195, -0.5029860367700724);
	arr[i++] = _newGrad4(0.4004672082940195, 0.7504883828755602, -0.15296486218853164, -0.5029860367700724);
	arr[i++] = _newGrad4(0.3239847771997537, 0.6740059517812944, 0.3239847771997537, -0.5794684678643381);
	arr[i++] = _newGrad4(0.5794684678643381, -0.6740059517812944, -0.3239847771997537, -0.3239847771997537);
	arr[i++] = _newGrad4(0.5029860367700724, -0.7504883828755602, -0.4004672082940195, 0.15296486218853164);
	arr[i++] = _newGrad4(0.5029860367700724, -0.7504883828755602, 0.15296486218853164, -0.4004672082940195);
	arr[i++] = _newGrad4(0.4553054119602712, -0.8828161875373585, 0.08164729285680945, 0.08164729285680945);
	arr[i++] = _newGrad4(0.8828161875373585, -0.4553054119602712, -0.08164729285680945, -0.08164729285680945);
	arr[i++] = _newGrad4(0.7504883828755602, -0.5029860367700724, -0.15296486218853164, 0.4004672082940195);
	arr[i++] = _newGrad4(0.7504883828755602, -0.5029860367700724, 0.4004672082940195, -0.15296486218853164);
	arr[i++] = _newGrad4(0.6740059517812944, -0.5794684678643381, 0.3239847771997537, 0.3239847771997537);
	arr[i++] = _newGrad4(0.5794684678643381, -0.3239847771997537, -0.6740059517812944, -0.3239847771997537);
	arr[i++] = _newGrad4(0.5029860367700724, -0.4004672082940195, -0.7504883828755602, 0.15296486218853164);
	arr[i++] = _newGrad4(0.5029860367700724, 0.15296486218853164, -0.7504883828755602, -0.4004672082940195);
	arr[i++] = _newGrad4(0.4553054119602712, 0.08164729285680945, -0.8828161875373585, 0.08164729285680945);
	arr[i++] = _newGrad4(0.8828161875373585, -0.08164729285680945, -0.4553054119602712, -0.08164729285680945);
	arr[i++] = _newGrad4(0.7504883828755602, -0.15296486218853164, -0.5029860367700724, 0.4004672082940195);
	arr[i++] = _newGrad4(0.7504883828755602, 0.4004672082940195, -0.5029860367700724, -0.15296486218853164);
	arr[i++] = _newGrad4(0.6740059517812944, 0.3239847771997537, -0.5794684678643381, 0.3239847771997537);
	arr[i++] = _newGrad4(0.5794684678643381, -0.3239847771997537, -0.3239847771997537, -0.6740059517812944);
	arr[i++] = _newGrad4(0.5029860367700724, -0.4004672082940195, 0.15296486218853164, -0.7504883828755602);
	arr[i++] = _newGrad4(0.5029860367700724, 0.15296486218853164, -0.4004672082940195, -0.7504883828755602);
	arr[i++] = _newGrad4(0.4553054119602712, 0.08164729285680945, 0.08164729285680945, -0.8828161875373585);
	arr[i++] = _newGrad4(0.8828161875373585, -0.08164729285680945, -0.08164729285680945, -0.4553054119602712);
	arr[i++] = _newGrad4(0.7504883828755602, -0.15296486218853164, 0.4004672082940195, -0.5029860367700724);
	arr[i++] = _newGrad4(0.7504883828755602, 0.4004672082940195, -0.15296486218853164, -0.5029860367700724);
	arr[i++] = _newGrad4(0.6740059517812944, 0.3239847771997537, 0.3239847771997537, -0.5794684678643381);
	arr[i++] = _newGrad4(0.03381941603233842, 0.03381941603233842, 0.03381941603233842, 0.9982828964265062);
	arr[i++] = _newGrad4(-0.044802370851755174, -0.044802370851755174, 0.508629699630796, 0.8586508742123365);
	arr[i++] = _newGrad4(-0.044802370851755174, 0.508629699630796, -0.044802370851755174, 0.8586508742123365);
	arr[i++] = _newGrad4(-0.12128480194602098, 0.4321472685365301, 0.4321472685365301, 0.7821684431180708);
	arr[i++] = _newGrad4(0.508629699630796, -0.044802370851755174, -0.044802370851755174, 0.8586508742123365);
	arr[i++] = _newGrad4(0.4321472685365301, -0.12128480194602098, 0.4321472685365301, 0.7821684431180708);
	arr[i++] = _newGrad4(0.4321472685365301, 0.4321472685365301, -0.12128480194602098, 0.7821684431180708);
	arr[i++] = _newGrad4(0.37968289875261624, 0.37968289875261624, 0.37968289875261624, 0.753341017856078);
	arr[i++] = _newGrad4(0.03381941603233842, 0.03381941603233842, 0.9982828964265062, 0.03381941603233842);
	arr[i++] = _newGrad4(-0.044802370851755174, 0.044802370851755174, 0.8586508742123365, 0.508629699630796);
	arr[i++] = _newGrad4(-0.044802370851755174, 0.508629699630796, 0.8586508742123365, -0.044802370851755174);
	arr[i++] = _newGrad4(-0.12128480194602098, 0.4321472685365301, 0.7821684431180708, 0.4321472685365301);
	arr[i++] = _newGrad4(0.508629699630796, -0.044802370851755174, 0.8586508742123365, -0.044802370851755174);
	arr[i++] = _newGrad4(0.4321472685365301, -0.12128480194602098, 0.7821684431180708, 0.4321472685365301);
	arr[i++] = _newGrad4(0.4321472685365301, 0.4321472685365301, 0.7821684431180708, -0.12128480194602098);
	arr[i++] = _newGrad4(0.37968289875261624, 0.37968289875261624, 0.753341017856078, 0.37968289875261624);
	arr[i++] = _newGrad4(0.03381941603233842, 0.9982828964265062, 0.03381941603233842, 0.03381941603233842);
	arr[i++] = _newGrad4(-0.044802370851755174, 0.8586508742123365, -0.044802370851755174, 0.508629699630796);
	arr[i++] = _newGrad4(-0.044802370851755174, 0.8586508742123365, 0.508629699630796, -0.044802370851755174);
	arr[i++] = _newGrad4(-0.12128480194602098, 0.7821684431180708, 0.4321472685365301, 0.4321472685365301);
	arr[i++] = _newGrad4(0.508629699630796, 0.8586508742123365, -0.044802370851755174, -0.044802370851755174);
	arr[i++] = _newGrad4(0.4321472685365301, 0.7821684431180708, -0.12128480194602098, 0.4321472685365301);
	arr[i++] = _newGrad4(0.4321472685365301, 0.7821684431180708, 0.4321472685365301, -0.12128480194602098);
	arr[i++] = _newGrad4(0.37968289875261624, 0.753341017856078, 0.37968289875261624, 0.37968289875261624);
	arr[i++] = _newGrad4(0.9982828964265062, 0.03381941603233842, 0.03381941603233842, 0.03381941603233842);
	arr[i++] = _newGrad4(0.8586508742123365, -0.044802370851755174, -0.044802370851755174, 0.508629699630796);
	arr[i++] = _newGrad4(0.8586508742123365, -0.044802370851755174, 0.508629699630796, -0.044802370851755174);
	arr[i++] = _newGrad4(0.7821684431180708, -0.12128480194602098, 0.4321472685365301, 0.4321472685365301);
	arr[i++] = _newGrad4(0.8586508742123365, 0.508629699630796, -0.044802370851755174, -0.044802370851755174);
	arr[i++] = _newGrad4(0.7821684431180708, 0.4321472685365301, -0.12128480194602098, 0.4321472685365301);
	arr[i++] = _newGrad4(0.7821684431180708, 0.4321472685365301, 0.4321472685365301, -0.12128480194602098);
	arr[i++] = _newGrad4(0.753341017856078, 0.37968289875261624, 0.37968289875261624, 0.37968289875261624);
	Grad4 *gradients4D = _newGrad4Arr(PSIZE);
	for (int i = 0; i < 160; i++){
		arr[i].dx /= N4;
		arr[i].dy /= N4;
		arr[i].dz /= N4;
		arr[i].dw /= N4;
	}
	for (int i = 0; i < PSIZE; i++){
		gradients4D[i] = arr[i % 160];
	}
	return gradients4D;
}

OpenSimplexEnv* initOpenSimplex(){
	OpenSimplexEnv *ose = (OpenSimplexEnv *) malloc(sizeof(OpenSimplexEnv));
	ose->GRADIENTS_2D = _newGrad2ConstArray();
	ose->GRADIENTS_3D = _newGrad3ConstArray();
	ose->GRADIENTS_4D = _newGrad4ConstArray();
	ose->LOOKUP_2D = _newLatticePoint2DConstArray();
	ose->LOOKUP_3D = _newLatticePoint3DConstArray();
	ose->LOOKUP_4D = _newLatticePoint4DConstArray();
	return ose;
}

OpenSimplexGradients* newOpenSimplexGradients(OpenSimplexEnv *ose, long seed){
	OpenSimplexGradients *osg = (OpenSimplexGradients *)malloc(sizeof(OpenSimplexGradients));
	osg->perm = _newShortArr(PSIZE);
	osg->permGrad2 = _newGrad2Arr(PSIZE);
	osg->permGrad3 = _newGrad3Arr(PSIZE);
	osg->permGrad4 = _newGrad4Arr(PSIZE);
	short *source = _newShortArr(PSIZE);
	for (short i = 0; i < PSIZE; i++)
		source[i] = i;
	for (int i = PSIZE - 1; i >= 0; i--){
		seed = seed * 6364136223846793005L + 1442695040888963407L;
		int r = (int)((seed + 31) % (i + 1));
		if (r < 0)
			r += (i + 1);
		osg->perm[i] = source[r];
		osg->permGrad2[i] = ose->GRADIENTS_2D[osg->perm[i]];
		osg->permGrad3[i] = ose->GRADIENTS_3D[osg->perm[i]];
		osg->permGrad4[i] = ose->GRADIENTS_4D[osg->perm[i]];
		source[r] = source[i];
	}
	return osg;
}