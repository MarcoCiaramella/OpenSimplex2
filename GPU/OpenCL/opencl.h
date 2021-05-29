#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>




typedef struct {
    cl_device_id device;
    cl_context context;
    cl_program program;
    cl_command_queue queue;
    unsigned int width;
	unsigned int height;
} OpenCLEnv;

OpenSimplexEnv *initOpenSimplex();
OpenSimplexGradients *newOpenSimplexGradients(OpenSimplexEnv *ose, cl_long seed);
OpenCLEnv initOpenCL(char *kernel_filename, unsigned int width, unsigned int height);
void releaseOpenCL(OpenCLEnv* openCLEnv);