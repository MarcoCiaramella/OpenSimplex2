#include "bitmap.h"
#include "OpenSimplex2S.h"
//#include "OpenSimplex2F.h"



#define WIDTH 4096
#define HEIGHT 4096




int main(){
     OpenSimplexEnv* ose = initOpenSimplex();
     OpenSimplexGradients* osg = newOpenSimplexGradients(ose, 1234);
     //OpenCLEnv openCLEnv = initOpenCL("OpenSimplex2F\\OpenSimplex2F.cl", WIDTH, HEIGHT);
     OpenCLEnv openCLEnv = initOpenCL("OpenSimplex2S\\OpenSimplex2S.cl", WIDTH, HEIGHT);

     save_bitmap("img/noise2.bmp", WIDTH, HEIGHT, noise2(&openCLEnv, ose, osg));
     save_bitmap("img/noise2_XBeforeY.bmp", WIDTH, HEIGHT, noise2_XBeforeY(&openCLEnv, ose, osg));
     save_bitmap("img/noise3_Classic.bmp", WIDTH, HEIGHT, noise3_Classic(&openCLEnv, ose, osg));
     save_bitmap("img/noise3_XYBeforeZ.bmp", WIDTH, HEIGHT, noise3_XYBeforeZ(&openCLEnv, ose, osg));
     save_bitmap("img/noise3_XZBeforeY.bmp", WIDTH, HEIGHT, noise3_XZBeforeY(&openCLEnv, ose, osg));
     save_bitmap("img/noise4_Classic.bmp", WIDTH, HEIGHT, noise4_Classic(&openCLEnv, ose, osg));
     save_bitmap("img/noise4_XYBeforeZW.bmp", WIDTH, HEIGHT, noise4_XYBeforeZW(&openCLEnv, ose, osg));
     save_bitmap("img/noise4_XZBeforeYW.bmp", WIDTH, HEIGHT, noise4_XZBeforeYW(&openCLEnv, ose, osg));
     save_bitmap("img/noise4_XYZBeforeW.bmp", WIDTH, HEIGHT, noise4_XYZBeforeW(&openCLEnv, ose, osg));
     return 0;
}