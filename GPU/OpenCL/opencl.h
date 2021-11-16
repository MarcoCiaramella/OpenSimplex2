#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>




typedef struct {
    cl_device_id device;
    cl_context context;
    cl_program program;
    cl_command_queue queue;
} OpenCLEnv;

OpenCLEnv loadOpenCL(char *kernel_filename, unsigned int width, unsigned int height);
void releaseOpenCL(OpenCLEnv* openCLEnv);
double *run_kernel(
	OpenCLEnv* openCLEnv,
	char *function,
	void *host_ptr1,
	void *host_ptr2,
	void *host_ptr3,
	size_t size1,
	size_t size2,
	size_t size3);