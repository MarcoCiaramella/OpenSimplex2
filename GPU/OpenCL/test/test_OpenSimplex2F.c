#include "bitmap.h"
#include "../OpenSimplex2F.h"



#define SIZE_2D 4096
#define SIZE_3D 256
#define SIZE_4D 64
#define OFF_X 2048
#define OFF_Y 2048
#define OFF_Z 2048
#define OFF_W 2048
#define PERIOD 64.0
#define FREQ 1.0 / PERIOD


double* new_grid2D(size_t* size){
     size_t num_points = SIZE_2D*SIZE_2D;
     *size = num_points * 2 * sizeof(double);
     double* grid = (double*) malloc(*size);
     int i = 0;
     for (int y = 0; y < SIZE_2D; y++){
          for (int x = 0; x < SIZE_2D; x++){
               grid[i++] = (x + OFF_X) * FREQ;
               grid[i++] = (y + OFF_Y) * FREQ;
          }
     }
     return grid;
}

double* new_grid3D(size_t* size){
     size_t num_points = SIZE_3D*SIZE_3D*SIZE_3D;
     *size = num_points * 3 * sizeof(double);
     double* grid = (double*) malloc(*size);
     int i = 0;
     for (int z = 0; z < SIZE_3D; z++){
          for (int y = 0; y < SIZE_3D; y++){
               for (int x = 0; x < SIZE_3D; x++){
                    grid[i++] = (x + OFF_X) * FREQ;
                    grid[i++] = (y + OFF_Y) * FREQ;
                    grid[i++] = (z + OFF_Z) * FREQ;
               }
          }
     }
     return grid;
}
double* new_grid4D(size_t* size){
     size_t num_points = SIZE_4D*SIZE_4D*SIZE_4D*SIZE_4D;
     *size = num_points * 4 * sizeof(double);
     double* grid = (double*) malloc(*size);
     int i = 0;
     for (int w = 0; w < SIZE_4D; w++){
          for (int z = 0; z < SIZE_4D; z++){
               for (int y = 0; y < SIZE_4D; y++){
                    for (int x = 0; x < SIZE_4D; x++){
                         grid[i++] = (x + OFF_X) * FREQ;
                         grid[i++] = (y + OFF_Y) * FREQ;
                         grid[i++] = (z + OFF_Z) * FREQ;
                         grid[i++] = (w + OFF_W) * FREQ;
                    }
               }
          }
     }
     return grid;
}



int main(){
     OpenSimplexEnv* ose = initOpenSimplex();
     OpenSimplexGradients* osg = newOpenSimplexGradients(ose, 1234);
     OpenCLEnv openCLEnv = loadOpenCL("../OpenSimplex2F.cl");

     double* output_buffer;
     double* grid;
     size_t size_input_buffer;

     grid = new_grid2D(&size_input_buffer);
     output_buffer = noise2(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise2.bmp", output_buffer, SIZE_2D, SIZE_2D);
     output_buffer = noise2_XBeforeY(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise2_XBeforeY.bmp", output_buffer, SIZE_2D, SIZE_2D);

     grid = new_grid3D(&size_input_buffer);
     output_buffer = noise3_Classic(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise3_Classic.bmp", output_buffer, SIZE_3D, SIZE_3D);
     output_buffer = noise3_XYBeforeZ(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise3_XYBeforeZ.bmp", output_buffer, SIZE_3D, SIZE_3D);
     output_buffer = noise3_XZBeforeY(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise3_XZBeforeY.bmp", output_buffer, SIZE_3D, SIZE_3D);

     grid = new_grid4D(&size_input_buffer);
     output_buffer = noise4_Classic(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_Classic.bmp", output_buffer, SIZE_4D, SIZE_4D);
     output_buffer = noise4_XYBeforeZW(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_XYBeforeZW.bmp", output_buffer, SIZE_4D, SIZE_4D);
     output_buffer = noise4_XZBeforeYW(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_XZBeforeYW.bmp", output_buffer, SIZE_4D, SIZE_4D);
     output_buffer = noise4_XYZBeforeW(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_XYZBeforeW.bmp", output_buffer, SIZE_4D, SIZE_4D);

     releaseOpenCL(&openCLEnv);
     return 0;
}