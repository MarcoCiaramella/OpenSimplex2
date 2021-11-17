#include "bitmap.h"
#include "../OpenSimplex2S.h"



#define WIDTH 8
#define HEIGHT 8
#define SIZE WIDTH*HEIGHT


double* new_grid2D(size_t num_points, size_t* size){
     *size = num_points * 2 * sizeof(double);
     double* grid = (double*) malloc(*size);
     for (int i = 0; i < num_points; i++){
          grid[i] = 0;
          grid[i++] = 0;
     }
     return grid;
}


int main(){
     OpenSimplexEnv* ose = initOpenSimplex();
     OpenSimplexGradients* osg = newOpenSimplexGradients(ose, 1234);
     OpenCLEnv openCLEnv = loadOpenCL("../OpenSimplex2S.cl");

     double *output_buffer;
     size_t size_input_buffer;

     output_buffer = noise2(&openCLEnv, ose, osg, new_grid2D(SIZE, &size_input_buffer), size_input_buffer);
     save_bitmap("OpenSimplex2S_noise2.bmp", WIDTH, HEIGHT, output_buffer);
     output_buffer = noise2_XBeforeY(&openCLEnv, ose, osg, new_grid2D(SIZE, &size_input_buffer), size_input_buffer);
     save_bitmap("OpenSimplex2S_noise2_XBeforeY.bmp", WIDTH, HEIGHT, output_buffer);
     /*output_buffer = noise3_Classic(&openCLEnv, ose, osg);
     save_bitmap("OpenSimplex2S_noise3_Classic.bmp", WIDTH, HEIGHT, output_buffer);
     output_buffer = noise3_XYBeforeZ(&openCLEnv, ose, osg);
     save_bitmap("OpenSimplex2S_noise3_XYBeforeZ.bmp", WIDTH, HEIGHT, output_buffer);
     output_buffer = noise3_XZBeforeY(&openCLEnv, ose, osg);
     save_bitmap("OpenSimplex2S_noise3_XZBeforeY.bmp", WIDTH, HEIGHT, output_buffer);
     output_buffer = noise4_Classic(&openCLEnv, ose, osg);
     save_bitmap("OpenSimplex2S_noise4_Classic.bmp", WIDTH, HEIGHT, output_buffer);
     output_buffer = noise4_XYBeforeZW(&openCLEnv, ose, osg);
     save_bitmap("OpenSimplex2S_noise4_XYBeforeZW.bmp", WIDTH, HEIGHT, output_buffer);
     output_buffer = noise4_XZBeforeYW(&openCLEnv, ose, osg);
     save_bitmap("OpenSimplex2S_noise4_XZBeforeYW.bmp", WIDTH, HEIGHT, output_buffer);
     output_buffer = noise4_XYZBeforeW(&openCLEnv, ose, osg);
     save_bitmap("OpenSimplex2S_noise4_XYZBeforeW.bmp", WIDTH, HEIGHT, output_buffer);*/

     releaseOpenCL(&openCLEnv);
     return 0;
}