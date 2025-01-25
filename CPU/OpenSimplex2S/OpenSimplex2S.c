#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "OpenSimplex2S.h"

#define PSIZE 2048
#define PMASK 2047
#define N2 0.05481866495625118
#define N3 0.2781926117527186
#define N4 0.11127401889945551

/*
 * Utility
 */

int _fastFloor(double x)
{
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
}

Grad2 *_newGrad2Arr(unsigned int size)
{
	return (Grad2 *)malloc(sizeof(Grad2) * size);
}

Grad3 *_newGrad3Arr(unsigned int size)
{
	return (Grad3 *)malloc(sizeof(Grad3) * size);
}

Grad4 *_newGrad4Arr(unsigned int size)
{
	return (Grad4 *)malloc(sizeof(Grad4) * size);
}

short *_newShortArr(unsigned int size)
{
	return (short *)malloc(sizeof(short) * size);
}

LatticePoint2D *_newLatticePoint2D(int xsv, int ysv)
{
	LatticePoint2D *plp2D = (LatticePoint2D *)malloc(sizeof(LatticePoint2D));
	plp2D->xsv = xsv;
	plp2D->ysv = ysv;
	double ssv = (xsv + ysv) * -0.211324865405187;
	plp2D->dx = -xsv - ssv;
	plp2D->dy = -ysv - ssv;
	return plp2D;
}

LatticePoint3D *_newLatticePoint3D(int xrv, int yrv, int zrv, int lattice)
{
	LatticePoint3D *plp3D = (LatticePoint3D *)malloc(sizeof(LatticePoint3D));
	plp3D->dxr = -xrv + lattice * 0.5;
	plp3D->dyr = -yrv + lattice * 0.5;
	plp3D->dzr = -zrv + lattice * 0.5;
	plp3D->xrv = xrv + lattice * 1024;
	plp3D->yrv = yrv + lattice * 1024;
	plp3D->zrv = zrv + lattice * 1024;
	return plp3D;
}

LatticePoint4D *_newLatticePoint4D(int xsv, int ysv, int zsv, int wsv)
{
	LatticePoint4D *plp4D = (LatticePoint4D *)malloc(sizeof(LatticePoint4D));
	plp4D->xsv = xsv;
	plp4D->ysv = ysv;
	plp4D->zsv = zsv;
	plp4D->wsv = wsv;
	double ssv = (xsv + ysv + zsv + wsv) * -0.138196601125011;
	plp4D->dx = -xsv - ssv;
	plp4D->dy = -ysv - ssv;
	plp4D->dz = -zsv - ssv;
	plp4D->dw = -wsv - ssv;
	return plp4D;
}

LatticePoint2D **_newLatticePoint2DConstArray()
{
	LatticePoint2D **plp2DArr = (LatticePoint2D **)malloc(sizeof(LatticePoint2D *) * 8 * 4);
	for (int i = 0; i < 8; i++)
	{
		int i1, j1, i2, j2;
		if ((i & 1) == 0)
		{
			if ((i & 2) == 0)
			{
				i1 = -1;
				j1 = 0;
			}
			else
			{
				i1 = 1;
				j1 = 0;
			}
			if ((i & 4) == 0)
			{
				i2 = 0;
				j2 = -1;
			}
			else
			{
				i2 = 0;
				j2 = 1;
			}
		}
		else
		{
			if ((i & 2) != 0)
			{
				i1 = 2;
				j1 = 1;
			}
			else
			{
				i1 = 0;
				j1 = 1;
			}
			if ((i & 4) != 0)
			{
				i2 = 1;
				j2 = 2;
			}
			else
			{
				i2 = 1;
				j2 = 0;
			}
		}
		plp2DArr[i * 4 + 0] = _newLatticePoint2D(0, 0);
		plp2DArr[i * 4 + 1] = _newLatticePoint2D(1, 1);
		plp2DArr[i * 4 + 2] = _newLatticePoint2D(i1, j1);
		plp2DArr[i * 4 + 3] = _newLatticePoint2D(i2, j2);
	}
	return plp2DArr;
}

LatticePoint3D **_newLatticePoint3DConstArray()
{
	LatticePoint3D **plp3DArr = (LatticePoint3D **)malloc(sizeof(LatticePoint3D *) * 8);
	for (int i = 0; i < 8; i++)
	{
		int i1, j1, k1, i2, j2, k2;
		i1 = (i >> 0) & 1;
		j1 = (i >> 1) & 1;
		k1 = (i >> 2) & 1;
		i2 = i1 ^ 1;
		j2 = j1 ^ 1;
		k2 = k1 ^ 1;

		// The two points within this octant, one from each of the two cubic half-lattices.
		LatticePoint3D *c0 = _newLatticePoint3D(i1, j1, k1, 0);
		LatticePoint3D *c1 = _newLatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1);

		// (1, 0, 0) vs (0, 1, 1) away from octant.
		LatticePoint3D *c2 = _newLatticePoint3D(i1 ^ 1, j1, k1, 0);
		LatticePoint3D *c3 = _newLatticePoint3D(i1, j1 ^ 1, k1 ^ 1, 0);

		// (1, 0, 0) vs (0, 1, 1) away from octant, on second half-lattice.
		LatticePoint3D *c4 = _newLatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1);
		LatticePoint3D *c5 = _newLatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + (k2 ^ 1), 1);

		// (0, 1, 0) vs (1, 0, 1) away from octant.
		LatticePoint3D *c6 = _newLatticePoint3D(i1, j1 ^ 1, k1, 0);
		LatticePoint3D *c7 = _newLatticePoint3D(i1 ^ 1, j1, k1 ^ 1, 0);

		// (0, 1, 0) vs (1, 0, 1) away from octant, on second half-lattice.
		LatticePoint3D *c8 = _newLatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1);
		LatticePoint3D *c9 = _newLatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + (k2 ^ 1), 1);

		// (0, 0, 1) vs (1, 1, 0) away from octant.
		LatticePoint3D *cA = _newLatticePoint3D(i1, j1, k1 ^ 1, 0);
		LatticePoint3D *cB = _newLatticePoint3D(i1 ^ 1, j1 ^ 1, k1, 0);

		// (0, 0, 1) vs (1, 1, 0) away from octant, on second half-lattice.
		LatticePoint3D *cC = _newLatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1);
		LatticePoint3D *cD = _newLatticePoint3D(i1 + (i2 ^ 1), j1 + (j2 ^ 1), k1 + k2, 1);

		// First two points are guaranteed.
		c0->nextOnFailure = c0->nextOnSuccess = c1;
		c1->nextOnFailure = c1->nextOnSuccess = c2;

		// If c2 is in range, then we know c3 and c4 are not.
		c2->nextOnFailure = c3;
		c2->nextOnSuccess = c5;
		c3->nextOnFailure = c4;
		c3->nextOnSuccess = c4;

		// If c4 is in range, then we know c5 is not.
		c4->nextOnFailure = c5;
		c4->nextOnSuccess = c6;
		c5->nextOnFailure = c5->nextOnSuccess = c6;

		// If c6 is in range, then we know c7 and c8 are not.
		c6->nextOnFailure = c7;
		c6->nextOnSuccess = c9;
		c7->nextOnFailure = c8;
		c7->nextOnSuccess = c8;

		// If c8 is in range, then we know c9 is not.
		c8->nextOnFailure = c9;
		c8->nextOnSuccess = cA;
		c9->nextOnFailure = c9->nextOnSuccess = cA;

		// If cA is in range, then we know cB and cC are not.
		cA->nextOnFailure = cB;
		cA->nextOnSuccess = cD;
		cB->nextOnFailure = cC;
		cB->nextOnSuccess = cC;

		// If cC is in range, then we know cD is not.
		cC->nextOnFailure = cD;
		cC->nextOnSuccess = NULL;
		cD->nextOnFailure = cD->nextOnSuccess = NULL;

		plp3DArr[i] = c0;
	}
	return plp3DArr;
}

vect _newVect(unsigned int length, int *data)
{
	vect v = {length, data};
	return v;
}

vect *_newLatticePoint4DConstArray()
{
	vect *plp4DArr = (vect *)malloc(sizeof(vect) * 256);
	vect lookup4DPregen[256];
	int i = 0;
	int arr1[] = {0x15, 0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr1 / sizeof *arr1, arr1);

	int arr2[] = {0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr2 / sizeof *arr2, arr2);

	int arr3[] = {0x01, 0x05, 0x11, 0x15, 0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr3 / sizeof *arr3, arr3);

	int arr4[] = {0x01, 0x15, 0x16, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr4 / sizeof *arr4, arr4);

	int arr5[] = {0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr5 / sizeof *arr5, arr5);

	int arr6[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr6 / sizeof *arr6, arr6);

	int arr7[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr7 / sizeof *arr7, arr7);

	int arr8[] = {0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr8 / sizeof *arr8, arr8);

	int arr9[] = {0x04, 0x05, 0x14, 0x15, 0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr9 / sizeof *arr9, arr9);

	int arr10[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr10 / sizeof *arr10, arr10);

	int arr11[] = {0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr11 / sizeof *arr11, arr11);

	int arr12[] = {0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x6A, 0x9A, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr12 / sizeof *arr12, arr12);

	int arr13[] = {0x04, 0x15, 0x19, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr13 / sizeof *arr13, arr13);

	int arr14[] = {0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr14 / sizeof *arr14, arr14);

	int arr15[] = {0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x5E, 0x6A, 0x9A, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr15 / sizeof *arr15, arr15);

	int arr16[] = {0x05, 0x15, 0x1A, 0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x5E, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr16 / sizeof *arr16, arr16);

	int arr17[] = {0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr17 / sizeof *arr17, arr17);

	int arr18[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr18 / sizeof *arr18, arr18);

	int arr19[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr19 / sizeof *arr19, arr19);

	int arr20[] = {0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr20 / sizeof *arr20, arr20);

	int arr21[] = {0x14, 0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr21 / sizeof *arr21, arr21);

	int arr22[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr22 / sizeof *arr22, arr22);

	int arr23[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr23 / sizeof *arr23, arr23);

	int arr24[] = {0x15, 0x16, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x6B, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr24 / sizeof *arr24, arr24);

	int arr25[] = {0x14, 0x15, 0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr25 / sizeof *arr25, arr25);

	int arr26[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr26 / sizeof *arr26, arr26);

	int arr27[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr27 / sizeof *arr27, arr27);

	int arr28[] = {0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x6B, 0x9A, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr28 / sizeof *arr28, arr28);

	int arr29[] = {0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr29 / sizeof *arr29, arr29);

	int arr30[] = {0x15, 0x19, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr30 / sizeof *arr30, arr30);

	int arr31[] = {0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x9A, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr31 / sizeof *arr31, arr31);

	int arr32[] = {0x15, 0x1A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr32 / sizeof *arr32, arr32);

	int arr33[] = {0x10, 0x11, 0x14, 0x15, 0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr33 / sizeof *arr33, arr33);

	int arr34[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr34 / sizeof *arr34, arr34);

	int arr35[] = {0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr35 / sizeof *arr35, arr35);

	int arr36[] = {0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x67, 0x6A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr36 / sizeof *arr36, arr36);

	int arr37[] = {0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr37 / sizeof *arr37, arr37);

	int arr38[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr38 / sizeof *arr38, arr38);

	int arr39[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr39 / sizeof *arr39, arr39);

	int arr40[] = {0x15, 0x16, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x6B, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr40 / sizeof *arr40, arr40);

	int arr41[] = {0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr41 / sizeof *arr41, arr41);

	int arr42[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr42 / sizeof *arr42, arr42);

	int arr43[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr43 / sizeof *arr43, arr43);

	int arr44[] = {0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr44 / sizeof *arr44, arr44);

	int arr45[] = {0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x6D, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr45 / sizeof *arr45, arr45);

	int arr46[] = {0x15, 0x19, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x6E, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr46 / sizeof *arr46, arr46);

	int arr47[] = {0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr47 / sizeof *arr47, arr47);

	int arr48[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr48 / sizeof *arr48, arr48);

	int arr49[] = {0x10, 0x15, 0x25, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr49 / sizeof *arr49, arr49);

	int arr50[] = {0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr50 / sizeof *arr50, arr50);

	int arr51[] = {0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0x76, 0xA6, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr51 / sizeof *arr51, arr51);

	int arr52[] = {0x11, 0x15, 0x26, 0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x67, 0x6A, 0x76, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr52 / sizeof *arr52, arr52);

	int arr53[] = {0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr53 / sizeof *arr53, arr53);

	int arr54[] = {0x15, 0x25, 0x55, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr54 / sizeof *arr54, arr54);

	int arr55[] = {0x15, 0x25, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA6, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr55 / sizeof *arr55, arr55);

	int arr56[] = {0x15, 0x26, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr56 / sizeof *arr56, arr56);

	int arr57[] = {0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0x79, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr57 / sizeof *arr57, arr57);

	int arr58[] = {0x15, 0x25, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr58 / sizeof *arr58, arr58);

	int arr59[] = {0x15, 0x25, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr59 / sizeof *arr59, arr59);

	int arr60[] = {0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr60 / sizeof *arr60, arr60);

	int arr61[] = {0x14, 0x15, 0x29, 0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0x6D, 0x79, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr61 / sizeof *arr61, arr61);

	int arr62[] = {0x15, 0x29, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr62 / sizeof *arr62, arr62);

	int arr63[] = {0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr63 / sizeof *arr63, arr63);

	int arr64[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x7A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF};
	lookup4DPregen[i++] = _newVect(sizeof arr64 / sizeof *arr64, arr64);

	int arr65[] = {0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr65 / sizeof *arr65, arr65);

	int arr66[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr66 / sizeof *arr66, arr66);

	int arr67[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr67 / sizeof *arr67, arr67);

	int arr68[] = {0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr68 / sizeof *arr68, arr68);

	int arr69[] = {0x44, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr69 / sizeof *arr69, arr69);

	int arr70[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr70 / sizeof *arr70, arr70);

	int arr71[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr71 / sizeof *arr71, arr71);

	int arr72[] = {0x45, 0x46, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr72 / sizeof *arr72, arr72);

	int arr73[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr73 / sizeof *arr73, arr73);

	int arr74[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr74 / sizeof *arr74, arr74);

	int arr75[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr75 / sizeof *arr75, arr75);

	int arr76[] = {0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x9A, 0x9B, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr76 / sizeof *arr76, arr76);

	int arr77[] = {0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr77 / sizeof *arr77, arr77);

	int arr78[] = {0x45, 0x49, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr78 / sizeof *arr78, arr78);

	int arr79[] = {0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x99, 0x9A, 0x9E, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr79 / sizeof *arr79, arr79);

	int arr80[] = {0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr80 / sizeof *arr80, arr80);

	int arr81[] = {0x50, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr81 / sizeof *arr81, arr81);

	int arr82[] = {0x51, 0x55, 0x56, 0x59, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr82 / sizeof *arr82, arr82);

	int arr83[] = {0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr83 / sizeof *arr83, arr83);

	int arr84[] = {0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr84 / sizeof *arr84, arr84);

	int arr85[] = {0x54, 0x55, 0x56, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr85 / sizeof *arr85, arr85);

	int arr86[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr86 / sizeof *arr86, arr86);

	int arr87[] = {0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr87 / sizeof *arr87, arr87);

	int arr88[] = {0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr88 / sizeof *arr88, arr88);

	int arr89[] = {0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr89 / sizeof *arr89, arr89);

	int arr90[] = {0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr90 / sizeof *arr90, arr90);

	int arr91[] = {0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr91 / sizeof *arr91, arr91);

	int arr92[] = {0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr92 / sizeof *arr92, arr92);

	int arr93[] = {0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr93 / sizeof *arr93, arr93);

	int arr94[] = {0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr94 / sizeof *arr94, arr94);

	int arr95[] = {0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr95 / sizeof *arr95, arr95);

	int arr96[] = {0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr96 / sizeof *arr96, arr96);

	int arr97[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr97 / sizeof *arr97, arr97);

	int arr98[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr98 / sizeof *arr98, arr98);

	int arr99[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr99 / sizeof *arr99, arr99);

	int arr100[] = {0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xA7, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr100 / sizeof *arr100, arr100);

	int arr101[] = {0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr101 / sizeof *arr101, arr101);

	int arr102[] = {0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr102 / sizeof *arr102, arr102);

	int arr103[] = {0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr103 / sizeof *arr103, arr103);

	int arr104[] = {0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr104 / sizeof *arr104, arr104);

	int arr105[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr105 / sizeof *arr105, arr105);

	int arr106[] = {0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr106 / sizeof *arr106, arr106);

	int arr107[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr107 / sizeof *arr107, arr107);

	int arr108[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr108 / sizeof *arr108, arr108);

	int arr109[] = {0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAD, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr109 / sizeof *arr109, arr109);

	int arr110[] = {0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr110 / sizeof *arr110, arr110);

	int arr111[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr111 / sizeof *arr111, arr111);

	int arr112[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr112 / sizeof *arr112, arr112);

	int arr113[] = {0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr113 / sizeof *arr113, arr113);

	int arr114[] = {0x51, 0x55, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr114 / sizeof *arr114, arr114);

	int arr115[] = {0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr115 / sizeof *arr115, arr115);

	int arr116[] = {0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x6A, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr116 / sizeof *arr116, arr116);

	int arr117[] = {0x54, 0x55, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr117 / sizeof *arr117, arr117);

	int arr118[] = {0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr118 / sizeof *arr118, arr118);

	int arr119[] = {0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr119 / sizeof *arr119, arr119);

	int arr120[] = {0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr120 / sizeof *arr120, arr120);

	int arr121[] = {0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr121 / sizeof *arr121, arr121);

	int arr122[] = {0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr122 / sizeof *arr122, arr122);

	int arr123[] = {0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr123 / sizeof *arr123, arr123);

	int arr124[] = {0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr124 / sizeof *arr124, arr124);

	int arr125[] = {0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr125 / sizeof *arr125, arr125);

	int arr126[] = {0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr126 / sizeof *arr126, arr126);

	int arr127[] = {0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr127 / sizeof *arr127, arr127);

	int arr128[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF};
	lookup4DPregen[i++] = _newVect(sizeof arr128 / sizeof *arr128, arr128);

	int arr129[] = {0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr129 / sizeof *arr129, arr129);

	int arr130[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr130 / sizeof *arr130, arr130);

	int arr131[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr131 / sizeof *arr131, arr131);

	int arr132[] = {0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr132 / sizeof *arr132, arr132);

	int arr133[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr133 / sizeof *arr133, arr133);

	int arr134[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr134 / sizeof *arr134, arr134);

	int arr135[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr135 / sizeof *arr135, arr135);

	int arr136[] = {0x45, 0x46, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr136 / sizeof *arr136, arr136);

	int arr137[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr137 / sizeof *arr137, arr137);

	int arr138[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr138 / sizeof *arr138, arr138);

	int arr139[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr139 / sizeof *arr139, arr139);

	int arr140[] = {0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr140 / sizeof *arr140, arr140);

	int arr141[] = {0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x95, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr141 / sizeof *arr141, arr141);

	int arr142[] = {0x45, 0x49, 0x55, 0x59, 0x5A, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr142 / sizeof *arr142, arr142);

	int arr143[] = {0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr143 / sizeof *arr143, arr143);

	int arr144[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr144 / sizeof *arr144, arr144);

	int arr145[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr145 / sizeof *arr145, arr145);

	int arr146[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr146 / sizeof *arr146, arr146);

	int arr147[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr147 / sizeof *arr147, arr147);

	int arr148[] = {0x51, 0x52, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr148 / sizeof *arr148, arr148);

	int arr149[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr149 / sizeof *arr149, arr149);

	int arr150[] = {0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr150 / sizeof *arr150, arr150);

	int arr151[] = {0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr151 / sizeof *arr151, arr151);

	int arr152[] = {0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr152 / sizeof *arr152, arr152);

	int arr153[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr153 / sizeof *arr153, arr153);

	int arr154[] = {0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr154 / sizeof *arr154, arr154);

	int arr155[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr155 / sizeof *arr155, arr155);

	int arr156[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr156 / sizeof *arr156, arr156);

	int arr157[] = {0x54, 0x55, 0x58, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr157 / sizeof *arr157, arr157);

	int arr158[] = {0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr158 / sizeof *arr158, arr158);

	int arr159[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr159 / sizeof *arr159, arr159);

	int arr160[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr160 / sizeof *arr160, arr160);

	int arr161[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr161 / sizeof *arr161, arr161);

	int arr162[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr162 / sizeof *arr162, arr162);

	int arr163[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr163 / sizeof *arr163, arr163);

	int arr164[] = {0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr164 / sizeof *arr164, arr164);

	int arr165[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr165 / sizeof *arr165, arr165);

	int arr166[] = {0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr166 / sizeof *arr166, arr166);

	int arr167[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr167 / sizeof *arr167, arr167);

	int arr168[] = {0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr168 / sizeof *arr168, arr168);

	int arr169[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr169 / sizeof *arr169, arr169);

	int arr170[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr170 / sizeof *arr170, arr170);

	int arr171[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA};
	lookup4DPregen[i++] = _newVect(sizeof arr171 / sizeof *arr171, arr171);

	int arr172[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB};
	lookup4DPregen[i++] = _newVect(sizeof arr172 / sizeof *arr172, arr172);

	int arr173[] = {0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr173 / sizeof *arr173, arr173);

	int arr174[] = {0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr174 / sizeof *arr174, arr174);

	int arr175[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAE};
	lookup4DPregen[i++] = _newVect(sizeof arr175 / sizeof *arr175, arr175);

	int arr176[] = {0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xAF};
	lookup4DPregen[i++] = _newVect(sizeof arr176 / sizeof *arr176, arr176);

	int arr177[] = {0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr177 / sizeof *arr177, arr177);

	int arr178[] = {0x51, 0x55, 0x61, 0x65, 0x66, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr178 / sizeof *arr178, arr178);

	int arr179[] = {0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr179 / sizeof *arr179, arr179);

	int arr180[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr180 / sizeof *arr180, arr180);

	int arr181[] = {0x54, 0x55, 0x64, 0x65, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr181 / sizeof *arr181, arr181);

	int arr182[] = {0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr182 / sizeof *arr182, arr182);

	int arr183[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr183 / sizeof *arr183, arr183);

	int arr184[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr184 / sizeof *arr184, arr184);

	int arr185[] = {0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr185 / sizeof *arr185, arr185);

	int arr186[] = {0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr186 / sizeof *arr186, arr186);

	int arr187[] = {0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr187 / sizeof *arr187, arr187);

	int arr188[] = {0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xBB};
	lookup4DPregen[i++] = _newVect(sizeof arr188 / sizeof *arr188, arr188);

	int arr189[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr189 / sizeof *arr189, arr189);

	int arr190[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr190 / sizeof *arr190, arr190);

	int arr191[] = {0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE};
	lookup4DPregen[i++] = _newVect(sizeof arr191 / sizeof *arr191, arr191);

	int arr192[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA};
	lookup4DPregen[i++] = _newVect(sizeof arr192 / sizeof *arr192, arr192);

	int arr193[] = {0x40, 0x45, 0x51, 0x54, 0x55, 0x85, 0x91, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr193 / sizeof *arr193, arr193);

	int arr194[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr194 / sizeof *arr194, arr194);

	int arr195[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xD6, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr195 / sizeof *arr195, arr195);

	int arr196[] = {0x41, 0x45, 0x51, 0x55, 0x56, 0x86, 0x92, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB, 0xD6, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr196 / sizeof *arr196, arr196);

	int arr197[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr197 / sizeof *arr197, arr197);

	int arr198[] = {0x45, 0x55, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xDA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr198 / sizeof *arr198, arr198);

	int arr199[] = {0x45, 0x55, 0x56, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xDA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr199 / sizeof *arr199, arr199);

	int arr200[] = {0x45, 0x55, 0x56, 0x86, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr200 / sizeof *arr200, arr200);

	int arr201[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xD9, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr201 / sizeof *arr201, arr201);

	int arr202[] = {0x45, 0x55, 0x59, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xDA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr202 / sizeof *arr202, arr202);

	int arr203[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xDA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr203 / sizeof *arr203, arr203);

	int arr204[] = {0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr204 / sizeof *arr204, arr204);

	int arr205[] = {0x44, 0x45, 0x54, 0x55, 0x59, 0x89, 0x95, 0x98, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE, 0xD9, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr205 / sizeof *arr205, arr205);

	int arr206[] = {0x45, 0x55, 0x59, 0x89, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr206 / sizeof *arr206, arr206);

	int arr207[] = {0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr207 / sizeof *arr207, arr207);

	int arr208[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xDA, 0xEA, 0xEF};
	lookup4DPregen[i++] = _newVect(sizeof arr208 / sizeof *arr208, arr208);

	int arr209[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr209 / sizeof *arr209, arr209);

	int arr210[] = {0x51, 0x55, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr210 / sizeof *arr210, arr210);

	int arr211[] = {0x51, 0x55, 0x56, 0x91, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr211 / sizeof *arr211, arr211);

	int arr212[] = {0x51, 0x55, 0x56, 0x92, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr212 / sizeof *arr212, arr212);

	int arr213[] = {0x54, 0x55, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr213 / sizeof *arr213, arr213);

	int arr214[] = {0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr214 / sizeof *arr214, arr214);

	int arr215[] = {0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr215 / sizeof *arr215, arr215);

	int arr216[] = {0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr216 / sizeof *arr216, arr216);

	int arr217[] = {0x54, 0x55, 0x59, 0x94, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr217 / sizeof *arr217, arr217);

	int arr218[] = {0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr218 / sizeof *arr218, arr218);

	int arr219[] = {0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr219 / sizeof *arr219, arr219);

	int arr220[] = {0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr220 / sizeof *arr220, arr220);

	int arr221[] = {0x54, 0x55, 0x59, 0x95, 0x98, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr221 / sizeof *arr221, arr221);

	int arr222[] = {0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr222 / sizeof *arr222, arr222);

	int arr223[] = {0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr223 / sizeof *arr223, arr223);

	int arr224[] = {0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xEA, 0xEF};
	lookup4DPregen[i++] = _newVect(sizeof arr224 / sizeof *arr224, arr224);

	int arr225[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xE5, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr225 / sizeof *arr225, arr225);

	int arr226[] = {0x51, 0x55, 0x65, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr226 / sizeof *arr226, arr226);

	int arr227[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr227 / sizeof *arr227, arr227);

	int arr228[] = {0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr228 / sizeof *arr228, arr228);

	int arr229[] = {0x54, 0x55, 0x65, 0x94, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr229 / sizeof *arr229, arr229);

	int arr230[] = {0x55, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr230 / sizeof *arr230, arr230);

	int arr231[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr231 / sizeof *arr231, arr231);

	int arr232[] = {0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr232 / sizeof *arr232, arr232);

	int arr233[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x94, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr233 / sizeof *arr233, arr233);

	int arr234[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr234 / sizeof *arr234, arr234);

	int arr235[] = {0x55, 0x56, 0x59, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr235 / sizeof *arr235, arr235);

	int arr236[] = {0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA, 0xEB};
	lookup4DPregen[i++] = _newVect(sizeof arr236 / sizeof *arr236, arr236);

	int arr237[] = {0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr237 / sizeof *arr237, arr237);

	int arr238[] = {0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr238 / sizeof *arr238, arr238);

	int arr239[] = {0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE};
	lookup4DPregen[i++] = _newVect(sizeof arr239 / sizeof *arr239, arr239);

	int arr240[] = {0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr240 / sizeof *arr240, arr240);

	int arr241[] = {0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA1, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA, 0xE5, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr241 / sizeof *arr241, arr241);

	int arr242[] = {0x51, 0x55, 0x65, 0x95, 0xA1, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr242 / sizeof *arr242, arr242);

	int arr243[] = {0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr243 / sizeof *arr243, arr243);

	int arr244[] = {0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xE6, 0xEA, 0xFB};
	lookup4DPregen[i++] = _newVect(sizeof arr244 / sizeof *arr244, arr244);

	int arr245[] = {0x54, 0x55, 0x65, 0x95, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr245 / sizeof *arr245, arr245);

	int arr246[] = {0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr246 / sizeof *arr246, arr246);

	int arr247[] = {0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr247 / sizeof *arr247, arr247);

	int arr248[] = {0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xEA, 0xFB};
	lookup4DPregen[i++] = _newVect(sizeof arr248 / sizeof *arr248, arr248);

	int arr249[] = {0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr249 / sizeof *arr249, arr249);

	int arr250[] = {0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr250 / sizeof *arr250, arr250);

	int arr251[] = {0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA};
	lookup4DPregen[i++] = _newVect(sizeof arr251 / sizeof *arr251, arr251);

	int arr252[] = {0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr252 / sizeof *arr252, arr252);

	int arr253[] = {0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xE9, 0xEA, 0xFE};
	lookup4DPregen[i++] = _newVect(sizeof arr253 / sizeof *arr253, arr253);

	int arr254[] = {0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA, 0xFE};
	lookup4DPregen[i++] = _newVect(sizeof arr254 / sizeof *arr254, arr254);

	int arr255[] = {0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr255 / sizeof *arr255, arr255);

	int arr256[] = {0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA, 0xEA};
	lookup4DPregen[i++] = _newVect(sizeof arr256 / sizeof *arr256, arr256);

	LatticePoint4D **latticePoints = (LatticePoint4D **)malloc(sizeof(LatticePoint4D *) * 256);
	for (int i = 0; i < 256; i++)
	{
		int cx = ((i >> 0) & 3) - 1;
		int cy = ((i >> 2) & 3) - 1;
		int cz = ((i >> 4) & 3) - 1;
		int cw = ((i >> 6) & 3) - 1;
		latticePoints[i] = _newLatticePoint4D(cx, cy, cz, cw);
	}
	for (int i = 0; i < 256; i++)
	{
		vect v;
		v.data = (LatticePoint4D *)malloc(sizeof(LatticePoint4D) * lookup4DPregen[i].length);
		v.length = lookup4DPregen[i].length;
		for (int j = 0; j < lookup4DPregen[i].length; j++)
		{
			((LatticePoint4D *)v.data)[j] = *(latticePoints[((int *)lookup4DPregen[i].data)[j]]);
		}
		plp4DArr[i] = v;
	}
	return plp4DArr;
}

Grad2 _newGrad2(double dx, double dy)
{
	Grad2 grad2;
	grad2.dx = dx;
	grad2.dy = dy;
	return grad2;
}

Grad3 _newGrad3(double dx, double dy, double dz)
{
	Grad3 grad3;
	grad3.dx = dx;
	grad3.dy = dy;
	grad3.dz = dz;
	return grad3;
}

Grad4 _newGrad4(double dx, double dy, double dz, double dw)
{
	Grad4 grad4;
	grad4.dx = dx;
	grad4.dy = dy;
	grad4.dz = dz;
	grad4.dw = dw;
	return grad4;
}

Grad2 *_newGrad2ConstArray()
{
	Grad2 *arr = (Grad2 *)malloc(sizeof(Grad2) * 24);
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
	for (int i = 0; i < 24; i++)
	{
		arr[i].dx /= N2;
		arr[i].dy /= N2;
	}
	for (int i = 0; i < PSIZE; i++)
	{
		gradients2D[i] = arr[i % 24];
	}
	return gradients2D;
}

Grad3 *_newGrad3ConstArray()
{
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
	for (int i = 0; i < 48; i++)
	{
		arr[i].dx /= N3;
		arr[i].dy /= N3;
		arr[i].dz /= N3;
	}
	for (int i = 0; i < PSIZE; i++)
	{
		gradients3D[i] = arr[i % 48];
	}
	return gradients3D;
}

Grad4 *_newGrad4ConstArray()
{
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
	for (int i = 0; i < 160; i++)
	{
		arr[i].dx /= N4;
		arr[i].dy /= N4;
		arr[i].dz /= N4;
		arr[i].dw /= N4;
	}
	for (int i = 0; i < PSIZE; i++)
	{
		gradients4D[i] = arr[i % 160];
	}
	return gradients4D;
}

/*
 * Noise Evaluators
 */

/**
 * 2D SuperSimplex noise base.
 * Lookup table implementation inspired by DigitalShadow.
 */
double _noise2_Base(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xs, double ys)
{
	double value = 0;

	// Get base points and offsets
	int xsb = _fastFloor(xs), ysb = _fastFloor(ys);
	double xsi = xs - xsb, ysi = ys - ysb;

	// Index to point list
	int a = (int)(xsi + ysi);
	int index =
		(a << 2) |
		(int)(xsi - ysi / 2 + 1 - a / 2.0) << 3 |
		(int)(ysi - xsi / 2 + 1 - a / 2.0) << 4;

	double ssi = (xsi + ysi) * -0.211324865405187;
	double xi = xsi + ssi, yi = ysi + ssi;

	// Point contributions
	for (int i = 0; i < 4; i++)
	{
		LatticePoint2D *c = ose->LOOKUP_2D[index + i];

		double dx = xi + c->dx, dy = yi + c->dy;
		double attn = 2.0 / 3.0 - dx * dx - dy * dy;
		if (attn <= 0)
			continue;

		int pxm = (xsb + c->xsv) & PMASK, pym = (ysb + c->ysv) & PMASK;
		Grad2 grad = osg->permGrad2[osg->perm[pxm] ^ pym];
		double extrapolation = grad.dx * dx + grad.dy * dy;

		attn *= attn;
		value += attn * attn * extrapolation;
	}

	return value;
}

/**
 * 2D SuperSimplex noise, standard lattice orientation.
 */
double noise2(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y)
{

	// Get points for A2* lattice
	double s = 0.366025403784439 * (x + y);
	double xs = x + s, ys = y + s;

	return _noise2_Base(ose, osg, xs, ys);
}

/**
 * 2D SuperSimplex noise, with Y pointing down the main diagonal.
 * Might be better for a 2D sandbox style game, where Y is vertical.
 * Probably slightly less optimal for heightmaps or continent maps.
 */
double noise2_XBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y)
{

	// Skew transform and rotation baked into one.
	double xx = x * 0.7071067811865476;
	double yy = y * 1.224744871380249;

	return _noise2_Base(ose, osg, yy + xx, yy - xx);
}

/**
 * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
 * Lookup table implementation inspired by DigitalShadow.
 * It was actually faster to narrow down the points in the loop itself,
 * than to build up the index with enough info to isolate 8 points.
 */
double _noise3_BCC(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xr, double yr, double zr)
{

	// Get base and offsets inside cube of first lattice.
	int xrb = _fastFloor(xr), yrb = _fastFloor(yr), zrb = _fastFloor(zr);
	double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

	// Identify which octant of the cube we're in. This determines which cell
	// in the other cubic lattice we're in, and also narrows down one point on each.
	int xht = (int)(xri + 0.5), yht = (int)(yri + 0.5), zht = (int)(zri + 0.5);
	int index = (xht << 0) | (yht << 1) | (zht << 2);

	// Point contributions
	double value = 0;
	LatticePoint3D *c = ose->LOOKUP_3D[index];
	while (c != NULL)
	{
		double dxr = xri + c->dxr, dyr = yri + c->dyr, dzr = zri + c->dzr;
		double attn = 0.75 - dxr * dxr - dyr * dyr - dzr * dzr;
		if (attn < 0)
		{
			c = c->nextOnFailure;
		}
		else
		{
			int pxm = (xrb + c->xrv) & PMASK, pym = (yrb + c->yrv) & PMASK, pzm = (zrb + c->zrv) & PMASK;
			Grad3 grad = osg->permGrad3[osg->perm[osg->perm[pxm] ^ pym] ^ pzm];
			double extrapolation = grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

			attn *= attn;
			value += attn * attn * extrapolation;
			c = c->nextOnSuccess;
		}
	}
	return value;
}

/**
 * 3D Re-oriented 8-point BCC noise, classic orientation
 * Proper substitute for what 3D SuperSimplex would be,
 * in light of Forbidden Formulae.
 * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
 */
double noise3_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z)
{

	// Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
	// If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
	// Orthonormal rotation. Not a skew transform.
	double r = (2.0 / 3.0) * (x + y + z);
	double xr = r - x, yr = r - y, zr = r - z;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}

/**
 * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Y).
 * Recommended for 3D terrain and time-varied animations.
 * The Z coordinate should always be the "different" coordinate in your use case.
 * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
 * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
 * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
 */
double noise3_XYBeforeZ(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z)
{

	// Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
	// Orthonormal rotation. Not a skew transform.
	double xy = x + y;
	double s2 = xy * -0.211324865405187;
	double zz = z * 0.577350269189626;
	double xr = x + s2 - zz, yr = y + s2 - zz;
	double zr = xy * 0.577350269189626 + zz;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}

/**
 * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Z).
 * Recommended for 3D terrain and time-varied animations.
 * The Y coordinate should always be the "different" coordinate in your use case.
 * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
 * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
 * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
 */
double noise3_XZBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z)
{

	// Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
	// Orthonormal rotation. Not a skew transform.
	double xz = x + z;
	double s2 = xz * -0.211324865405187;
	double yy = y * 0.577350269189626;
	double xr = x + s2 - yy;
	double zr = z + s2 - yy;
	double yr = xz * 0.577350269189626 + yy;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}

/**
 * 4D SuperSimplex noise base.
 * Using ultra-simple 4x4x4x4 lookup partitioning.
 * This isn't as elegant or SIMD/GPU/etc. portable as other approaches,
 * but it does compete performance-wise with optimized OpenSimplex1.
 */
double _noise4_Base(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xs, double ys, double zs, double ws)
{
	double value = 0;

	// Get base points and offsets
	int xsb = _fastFloor(xs), ysb = _fastFloor(ys), zsb = _fastFloor(zs), wsb = _fastFloor(ws);
	double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

	// Unskewed offsets
	double ssi = (xsi + ysi + zsi + wsi) * -0.138196601125011;
	double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;

	int index = ((_fastFloor(xs * 4) & 3) << 0) | ((_fastFloor(ys * 4) & 3) << 2) | ((_fastFloor(zs * 4) & 3) << 4) | ((_fastFloor(ws * 4) & 3) << 6);

	// Point contributions
	LatticePoint4D *c = (LatticePoint4D *)ose->LOOKUP_4D[index].data;
	for (int i = 0; i < ose->LOOKUP_4D[index].length; i++)
	{
		double dx = xi + c->dx, dy = yi + c->dy, dz = zi + c->dz, dw = wi + c->dw;
		double attn = 0.8 - dx * dx - dy * dy - dz * dz - dw * dw;
		if (attn > 0)
		{
			attn *= attn;

			int pxm = (xsb + c->xsv) & PMASK, pym = (ysb + c->ysv) & PMASK;
			int pzm = (zsb + c->zsv) & PMASK, pwm = (wsb + c->wsv) & PMASK;
			Grad4 grad = osg->permGrad4[osg->perm[osg->perm[osg->perm[pxm] ^ pym] ^ pzm] ^ pwm];
			double extrapolation = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

			value += attn * attn * extrapolation;
		}
		c++;
	}
	return value;
}

/**
 * 4D SuperSimplex noise, classic lattice orientation.
 */
double noise4_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w)
{

	// Get points for A4 lattice
	double s = 0.309016994374947 * (x + y + z + w);
	double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}

/**
 * 4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
 * Recommended for noise(x, y, sin(time), cos(time)) trick.
 */
double noise4_XYBeforeZW(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w)
{

	double s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718;
	double t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942;
	double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}

/**
 * 4D SuperSimplex noise, with XZ and YW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
 */
double noise4_XZBeforeYW(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w)
{

	double s2 = (x + z) * -0.28522513987434876941 + (y + w) * 0.83897065470611435718;
	double t2 = (y + w) * 0.21939749883706435719 + (x + z) * -0.48214856493302476942;
	double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}

/**
 * 4D SuperSimplex noise, with XYZ oriented like noise3_Classic,
 * and W for an extra degree of freedom.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 */
double noise4_XYZBeforeW(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w)
{

	double xyz = x + y + z;
	double ww = w * 1.118033988749894;
	double s2 = xyz * -0.16666666666666666 + ww;
	double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}

OpenSimplexEnv *initOpenSimplex()
{
	OpenSimplexEnv *ose = (OpenSimplexEnv *)malloc(sizeof(OpenSimplexEnv));
	ose->GRADIENTS_2D = _newGrad2ConstArray();
	ose->GRADIENTS_3D = _newGrad3ConstArray();
	ose->GRADIENTS_4D = _newGrad4ConstArray();
	ose->LOOKUP_2D = _newLatticePoint2DConstArray();
	ose->LOOKUP_3D = _newLatticePoint3DConstArray();
	ose->LOOKUP_4D = _newLatticePoint4DConstArray();
	return ose;
}

OpenSimplexGradients *newOpenSimplexGradients(OpenSimplexEnv *ose, long seed)
{
	OpenSimplexGradients *osg = (OpenSimplexGradients *)malloc(sizeof(OpenSimplexGradients));
	osg->perm = _newShortArr(PSIZE);
	osg->permGrad2 = _newGrad2Arr(PSIZE);
	osg->permGrad3 = _newGrad3Arr(PSIZE);
	osg->permGrad4 = _newGrad4Arr(PSIZE);
	short *source = _newShortArr(PSIZE);
	for (short i = 0; i < PSIZE; i++)
		source[i] = i;
	for (int i = PSIZE - 1; i >= 0; i--)
	{
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