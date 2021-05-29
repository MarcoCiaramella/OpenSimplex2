#include "OpenSimplex2F.h"






Grad2 _newGrad2(cl_double dx, cl_double dy){
    Grad2 grad2;
    grad2.dx = dx;
    grad2.dy = dy;
    return grad2;
}

Grad3 _newGrad3(cl_double dx, cl_double dy, cl_double dz){
    Grad3 grad3;
    grad3.dx = dx;
    grad3.dy = dy;
    grad3.dz = dz;
    return grad3;
}

Grad4 _newGrad4(cl_double dx, cl_double dy, cl_double dz, cl_double dw){
    Grad4 grad4;
    grad4.dx = dx;
    grad4.dy = dy;
    grad4.dz = dz;
    grad4.dw = dw;
    return grad4;
}

void _loadGrad2ConstArray(OpenSimplexEnv *ose){
    Grad2 arr[24];
    cl_int i = 0;
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
    for (cl_int i = 0; i < 24; i++) {
        arr[i].dx /= N2;
        arr[i].dy /= N2;
    }
	for (cl_int i = 0; i < PSIZE; i++) {
        ose->GRADIENTS_2D[i] = arr[i % 24];
    }
}

void _loadGrad3ConstArray(OpenSimplexEnv *ose){
	Grad3 arr[48];
	cl_int i = 0;
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
	for (cl_int i = 0; i < 48; i++){
		arr[i].dx /= N3;
		arr[i].dy /= N3;
		arr[i].dz /= N3;
	}
	for (cl_int i = 0; i < PSIZE; i++){
		ose->GRADIENTS_3D[i] = arr[i % 48];
	}
}

void _loadGrad4ConstArray(OpenSimplexEnv *ose){
	Grad4 arr[160];
	cl_int i = 0;
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
	for (cl_int i = 0; i < 160; i++){
		arr[i].dx /= N4;
		arr[i].dy /= N4;
		arr[i].dz /= N4;
		arr[i].dw /= N4;
	}
	for (cl_int i = 0; i < PSIZE; i++){
		ose->GRADIENTS_4D[i] = arr[i % 160];
	}
}

LatticePoint2D _newLatticePoint2D(cl_int xsv, cl_int ysv){
	LatticePoint2D lp2D;
	lp2D.xsv = xsv;
	lp2D.ysv = ysv;
	cl_double ssv = (xsv + ysv) * -0.211324865405187;
	lp2D.dx = -xsv - ssv;
	lp2D.dy = -ysv - ssv;
	return lp2D;
}

LatticePoint3D _newLatticePoint3D(cl_int xrv, cl_int yrv, cl_int zrv, cl_int lattice){
	LatticePoint3D lp3D;
	lp3D.dxr = -xrv + lattice * 0.5;
	lp3D.dyr = -yrv + lattice * 0.5;
	lp3D.dzr = -zrv + lattice * 0.5;
	lp3D.xrv = xrv + lattice * 1024;
	lp3D.yrv = yrv + lattice * 1024;
	lp3D.zrv = zrv + lattice * 1024;
	lp3D.nextOnFailure = -1;
	lp3D.nextOnSuccess = -1;
	return lp3D;
}

LatticePoint4D _newLatticePoint4D(cl_int xsv, cl_int ysv, cl_int zsv, cl_int wsv){
	LatticePoint4D lp4D;
	lp4D.xsv = xsv + 409;
	lp4D.ysv = ysv + 409;
	lp4D.zsv = zsv + 409;
	lp4D.wsv = wsv + 409;
	cl_double ssv = (xsv + ysv + zsv + wsv) * 0.309016994374947;
	lp4D.dx = -xsv - ssv;
	lp4D.dy = -ysv - ssv;
	lp4D.dz = -zsv - ssv;
	lp4D.dw = -wsv - ssv;
	lp4D.xsi = 0.2 - xsv;
	lp4D.ysi = 0.2 - ysv;
	lp4D.zsi = 0.2 - zsv;
	lp4D.wsi = 0.2 - wsv;
	lp4D.ssiDelta = (0.8 - xsv - ysv - zsv - wsv) * 0.309016994374947;
	return lp4D;
}

void _loadLatticePoint2DConstArray(OpenSimplexEnv *ose){
	ose->LOOKUP_2D[0] = _newLatticePoint2D(1, 0);
	ose->LOOKUP_2D[1] = _newLatticePoint2D(0, 0);
	ose->LOOKUP_2D[2] = _newLatticePoint2D(1, 1);
	ose->LOOKUP_2D[3] = _newLatticePoint2D(0, 1);
}

void _loadLatticePoint3DConstArray(OpenSimplexEnv *ose){
	int j = 7;
	for (cl_int i = 0; i < 8; i++){
		cl_int i1, j1, k1, i2, j2, k2;
		i1 = (i >> 0) & 1;
		j1 = (i >> 1) & 1;
		k1 = (i >> 2) & 1;
		i2 = i1 ^ 1;
		j2 = j1 ^ 1;
		k2 = k1 ^ 1;

		// The two points within this octant, one from each of the two cubic half-lattices.
		LatticePoint3D c0 = _newLatticePoint3D(i1, j1, k1, 0);
		LatticePoint3D c1 = _newLatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1);

		// Each single step away on the first half-lattice.
		LatticePoint3D c2 = _newLatticePoint3D(i1 ^ 1, j1, k1, 0);
		LatticePoint3D c3 = _newLatticePoint3D(i1, j1 ^ 1, k1, 0);
		LatticePoint3D c4 = _newLatticePoint3D(i1, j1, k1 ^ 1, 0);

		// Each single step away on the second half-lattice.
		LatticePoint3D c5 = _newLatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1);
		LatticePoint3D c6 = _newLatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1);
		LatticePoint3D c7 = _newLatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1);

		// First two are guaranteed.
		c0.nextOnFailure = c0.nextOnSuccess = j+1;
		c1.nextOnFailure = c1.nextOnSuccess = j+2;

		// Once we find one on the first half-lattice, the rest are out.
		// In addition, knowing c2 rules out c5.
		c2.nextOnFailure = j+3;
		c2.nextOnSuccess = j+6;
		c3.nextOnFailure = j+4;
		c3.nextOnSuccess = j+5;
		c4.nextOnFailure = c4.nextOnSuccess = j+5;

		// Once we find one on the second half-lattice, the rest are out.
		c5.nextOnFailure = j+6;
		c5.nextOnSuccess = -1;
		c6.nextOnFailure = j+7;
		c6.nextOnSuccess = -1;
		c7.nextOnFailure = c7.nextOnSuccess = -1;

		ose->LOOKUP_3D[i] = c0;
		ose->LOOKUP_3D[++j] = c1;
		ose->LOOKUP_3D[++j] = c2;
		ose->LOOKUP_3D[++j] = c3;
		ose->LOOKUP_3D[++j] = c4;
		ose->LOOKUP_3D[++j] = c5;
		ose->LOOKUP_3D[++j] = c6;
		ose->LOOKUP_3D[++j] = c7;
	}
}

void _loadLatticePoint4DConstArray(OpenSimplexEnv *ose){
	for (cl_int i = 0; i < 16; i++) {
		ose->VERTICES_4D[i] = _newLatticePoint4D((i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1, (i >> 3) & 1);
	}
}

OpenSimplexEnv* initOpenSimplex(){
	OpenSimplexEnv* ose = (OpenSimplexEnv*) malloc(sizeof(OpenSimplexEnv));
	_loadGrad2ConstArray(ose);
	_loadGrad3ConstArray(ose);
	_loadGrad4ConstArray(ose);
	_loadLatticePoint2DConstArray(ose);
	_loadLatticePoint3DConstArray(ose);
	_loadLatticePoint4DConstArray(ose);
	return ose;
}

OpenSimplexGradients* newOpenSimplexGradients(OpenSimplexEnv *ose, cl_long seed){
    OpenSimplexGradients* osg = (OpenSimplexGradients*) malloc(sizeof(OpenSimplexGradients));
    cl_short source[PSIZE];
    for (cl_short i = 0; i < PSIZE; i++){
        source[i] = i;
    }
    for (cl_int i = PSIZE - 1; i >= 0; i--){
        seed = seed * 6364136223846793005L + 1442695040888963407L;
		cl_int r = (cl_int)((seed + 31) % (i + 1));
		if (r < 0){
            r += (i + 1);
        }
		osg->perm[i] = source[r];
		osg->permGrad2[i] = ose->GRADIENTS_2D[osg->perm[i]];
		osg->permGrad3[i] = ose->GRADIENTS_3D[osg->perm[i]];
		osg->permGrad4[i] = ose->GRADIENTS_4D[osg->perm[i]];
		source[r] = source[i];
    }
	return osg;
}

double *noise2(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise2",
		osg->perm,
		osg->permGrad2,
		ose->LOOKUP_2D,
		sizeof(osg->perm),
		sizeof(osg->permGrad2),
		sizeof(ose->LOOKUP_2D));
}

double *noise2_XBeforeY(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise2_XBeforeY",
		osg->perm,
		osg->permGrad2,
		ose->LOOKUP_2D,
		sizeof(osg->perm),
		sizeof(osg->permGrad2),
		sizeof(ose->LOOKUP_2D));
}

double *noise3_Classic(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise3_Classic",
		osg->perm,
		osg->permGrad3,
		ose->LOOKUP_3D,
		sizeof(osg->perm),
		sizeof(osg->permGrad3),
		sizeof(ose->LOOKUP_3D));
}

double *noise3_XYBeforeZ(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise3_XYBeforeZ",
		osg->perm,
		osg->permGrad3,
		ose->LOOKUP_3D,
		sizeof(osg->perm),
		sizeof(osg->permGrad3),
		sizeof(ose->LOOKUP_3D));
}

double *noise3_XZBeforeY(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise3_XZBeforeY",
		osg->perm,
		osg->permGrad3,
		ose->LOOKUP_3D,
		sizeof(osg->perm),
		sizeof(osg->permGrad3),
		sizeof(ose->LOOKUP_3D));
}

double *noise4_Classic(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise4_Classic",
		osg->perm,
		osg->permGrad4,
		ose->VERTICES_4D,
		sizeof(osg->perm),
		sizeof(osg->permGrad4),
		sizeof(ose->VERTICES_4D));
}

double *noise4_XYBeforeZW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise4_XYBeforeZW",
		osg->perm,
		osg->permGrad4,
		ose->VERTICES_4D,
		sizeof(osg->perm),
		sizeof(osg->permGrad4),
		sizeof(ose->VERTICES_4D));
}

double *noise4_XZBeforeYW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise4_XZBeforeYW",
		osg->perm,
		osg->permGrad4,
		ose->VERTICES_4D,
		sizeof(osg->perm),
		sizeof(osg->permGrad4),
		sizeof(ose->VERTICES_4D));
}

double *noise4_XYZBeforeW(OpenCLEnv* openCLEnv, OpenSimplexEnv *ose, OpenSimplexGradients *osg){
	return run_kernel(
		openCLEnv,
		"noise4_XYZBeforeW",
		osg->perm,
		osg->permGrad4,
		ose->VERTICES_4D,
		sizeof(osg->perm),
		sizeof(osg->permGrad4),
		sizeof(ose->VERTICES_4D));
}