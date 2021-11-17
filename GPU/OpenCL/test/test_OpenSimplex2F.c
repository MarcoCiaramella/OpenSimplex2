#include "bitmap.h"
#include "../OpenSimplex2F.h"



#define WIDTH 64
#define HEIGHT 64
#define DEPTH 64
#define TIME 64


double* new_grid2D(size_t* size){
     size_t num_points = WIDTH*HEIGHT;
     *size = num_points * 2 * sizeof(double);
     double* grid = (double*) malloc(*size);
     int i = 0;
     for (int y = 0; y < HEIGHT; y++){
          for (int x = 0; x < WIDTH; x++){
               grid[i++] = x;
               grid[i++] = y;
          }
     }
     return grid;
}

double* new_grid3D(size_t* size){
     size_t num_points = WIDTH*HEIGHT*DEPTH;
     *size = num_points * 3 * sizeof(double);
     double* grid = (double*) malloc(*size);
     int i = 0;
     for (int z = 0; z < DEPTH; z++){
          for (int y = 0; y < HEIGHT; y++){
               for (int x = 0; x < WIDTH; x++){
                    grid[i++] = x;
                    grid[i++] = y;
                    grid[i++] = z;
               }
          }
     }
     return grid;
}
double* new_grid4D(size_t* size){
     size_t num_points = WIDTH*HEIGHT*DEPTH*TIME;
     *size = num_points * 4 * sizeof(double);
     double* grid = (double*) malloc(*size);
     int i = 0;
     for (int w = 0; w < TIME; w++){
          for (int z = 0; z < DEPTH; z++){
               for (int y = 0; y < HEIGHT; y++){
                    for (int x = 0; x < WIDTH; x++){
                         grid[i++] = x;
                         grid[i++] = y;
                         grid[i++] = z;
                         grid[i++] = w;
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
     save_bitmap("OpenSimplex2F_noise2.bmp", output_buffer, WIDTH, HEIGHT);
     output_buffer = noise2_XBeforeY(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise2_XBeforeY.bmp", output_buffer, WIDTH, HEIGHT);

     grid = new_grid3D(&size_input_buffer);
     output_buffer = noise3_Classic(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise3_Classic.bmp", output_buffer, WIDTH, HEIGHT);
     output_buffer = noise3_XYBeforeZ(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise3_XYBeforeZ.bmp", output_buffer, WIDTH, HEIGHT);
     output_buffer = noise3_XZBeforeY(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise3_XZBeforeY.bmp", output_buffer, WIDTH, HEIGHT);

     grid = new_grid4D(&size_input_buffer);
     output_buffer = noise4_Classic(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_Classic.bmp", output_buffer, WIDTH, HEIGHT);
     output_buffer = noise4_XYBeforeZW(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_XYBeforeZW.bmp", output_buffer, WIDTH, HEIGHT);
     output_buffer = noise4_XZBeforeYW(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_XZBeforeYW.bmp", output_buffer, WIDTH, HEIGHT);
     output_buffer = noise4_XYZBeforeW(&openCLEnv, ose, osg, grid, size_input_buffer);
     save_bitmap("OpenSimplex2F_noise4_XYZBeforeW.bmp", output_buffer, WIDTH, HEIGHT);

     releaseOpenCL(&openCLEnv);
     return 0;
}