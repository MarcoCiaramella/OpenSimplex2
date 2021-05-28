#include <stdio.h>
#include "bitmap.h"
#include <sys/timeb.h>
#include "OpenSimplex2F.h"
#include <math.h>



#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


#define WIDTH 4096
#define HEIGHT 4096




float get_time_s(struct timeb start, struct timeb end){
    unsigned long long ms = (unsigned long long) (1000.0 * (end.time - start.time) + (end.millitm - start.millitm));
    return ms/1000.0;
}

void exit_on_error(cl_int res){
     if (res != CL_SUCCESS){
          printf("Error.\n");
          exit(-1);
     }
}

cl_uint get_num_GPU_devices(cl_platform_id platform){
     cl_uint num_devices;
     exit_on_error(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices));
     return num_devices;
}

void get_GPU_device(cl_platform_id platform, cl_device_id* device){
     if (get_num_GPU_devices(platform) > 0){
          cl_device_id devices[1];
          exit_on_error(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, devices, NULL));
          *device = devices[0];
          return;
     }
     printf("No GPU device found.\n");
     exit(0);
}

void get_GPU_platform(cl_platform_id* gpu_platform){
     cl_uint num_platforms;
     cl_platform_id* platforms;
     exit_on_error(clGetPlatformIDs(0, NULL, &num_platforms));
     platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * num_platforms);
     exit_on_error(clGetPlatformIDs(num_platforms, platforms, NULL));
     for (int i = 0; i < num_platforms; i++){
          char platform_name[255];
          size_t size;
          exit_on_error(clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(platform_name), platform_name, &size));
          printf(platform_name);
          printf("\n");
          if (strcmp(platform_name, "NVIDIA CUDA") == 0 || strcmp(platform_name, "AMD Accelerated Parallel Processing") == 0){
               *gpu_platform = platforms[i];
               return;
          }
    }
    printf("No GPU platform found.\n");
    exit(0);
}

cl_ulong get_GPU_mem(cl_device_id device){
     cl_ulong size;
     clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
     return size;
}

char* read_file(char* filename){
     FILE* infile;
     char* buffer;
     long numbytes;

     infile = fopen(filename, "r");

     if (infile == NULL){
          printf("Error reading file.\n");
          exit(-1);
     }

     fseek(infile, 0L, SEEK_END);
     numbytes = ftell(infile);

     fseek(infile, 0L, SEEK_SET);

     buffer = (char*)calloc(numbytes, sizeof(char));

     if (buffer == NULL){
          printf("Error reading file.\n");
          exit(-1);
     }

     fread(buffer, sizeof(char), numbytes, infile);
     fclose(infile);

     return buffer;
}

char *get_error(cl_int res){
     switch (res){
     case CL_SUCCESS:
          return "CL_SUCCESS";
     case CL_DEVICE_NOT_FOUND:
          return "CL_DEVICE_NOT_FOUND";
     case CL_DEVICE_NOT_AVAILABLE:
          return "CL_DEVICE_NOT_AVAILABLE";
     case CL_COMPILER_NOT_AVAILABLE:
          return "CL_COMPILER_NOT_AVAILABLE";
     case CL_MEM_OBJECT_ALLOCATION_FAILURE:
          return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
     case CL_OUT_OF_RESOURCES:
          return "CL_OUT_OF_RESOURCES";
     case CL_OUT_OF_HOST_MEMORY:
          return "CL_OUT_OF_HOST_MEMORY";
     case CL_PROFILING_INFO_NOT_AVAILABLE:
          return "CL_PROFILING_INFO_NOT_AVAILABLE";
     case CL_MEM_COPY_OVERLAP:
          return "CL_MEM_COPY_OVERLAP";
     case CL_IMAGE_FORMAT_MISMATCH:
          return "CL_IMAGE_FORMAT_MISMATCH";
     case CL_IMAGE_FORMAT_NOT_SUPPORTED:
          return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
     case CL_BUILD_PROGRAM_FAILURE:
          return "CL_BUILD_PROGRAM_FAILURE";
     case CL_MAP_FAILURE:
          return "CL_MAP_FAILURE";
     case CL_MISALIGNED_SUB_BUFFER_OFFSET:
          return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
     case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
          return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
     case CL_COMPILE_PROGRAM_FAILURE:
          return "CL_COMPILE_PROGRAM_FAILURE";
     case CL_LINKER_NOT_AVAILABLE:
          return "CL_LINKER_NOT_AVAILABLE";
     case CL_LINK_PROGRAM_FAILURE:
          return "CL_LINK_PROGRAM_FAILURE";
     case CL_DEVICE_PARTITION_FAILED:
          return "CL_DEVICE_PARTITION_FAILED";
     case CL_KERNEL_ARG_INFO_NOT_AVAILABLE:
          return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
     case CL_INVALID_VALUE:
          return "CL_INVALID_VALUE";
     case CL_INVALID_DEVICE_TYPE:
          return "CL_INVALID_DEVICE_TYPE";
     case CL_INVALID_PLATFORM:
          return "CL_INVALID_PLATFORM";
     case CL_INVALID_DEVICE:
          return "CL_INVALID_DEVICE";
     case CL_INVALID_CONTEXT:
          return "CL_INVALID_CONTEXT";
     case CL_INVALID_QUEUE_PROPERTIES:
          return "CL_INVALID_QUEUE_PROPERTIES";
     case CL_INVALID_COMMAND_QUEUE:
          return "CL_INVALID_COMMAND_QUEUE";
     case CL_INVALID_HOST_PTR:
          return "CL_INVALID_HOST_PTR";
     case CL_INVALID_MEM_OBJECT:
          return "CL_INVALID_MEM_OBJECT";
     case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
          return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
     case CL_INVALID_IMAGE_SIZE:
          return "CL_INVALID_IMAGE_SIZE";
     case CL_INVALID_SAMPLER:
          return "CL_INVALID_SAMPLER";
     case CL_INVALID_BINARY:
          return "CL_INVALID_BINARY";
     case CL_INVALID_BUILD_OPTIONS:
          return "CL_INVALID_BUILD_OPTIONS";
     case CL_INVALID_PROGRAM:
          return "CL_INVALID_PROGRAM";
     case CL_INVALID_PROGRAM_EXECUTABLE:
          return "CL_INVALID_PROGRAM_EXECUTABLE";
     case CL_INVALID_KERNEL_NAME:
          return "CL_INVALID_KERNEL_NAME";
     case CL_INVALID_KERNEL_DEFINITION:
          return "CL_INVALID_KERNEL_DEFINITION";
     case CL_INVALID_KERNEL:
          return "CL_INVALID_KERNEL";
     case CL_INVALID_ARG_INDEX:
          return "CL_INVALID_ARG_INDEX";
     case CL_INVALID_ARG_VALUE:
          return "CL_INVALID_ARG_VALUE";
     case CL_INVALID_ARG_SIZE:
          return "CL_INVALID_ARG_SIZE";
     case CL_INVALID_KERNEL_ARGS:
          return "CL_INVALID_KERNEL_ARGS";
     case CL_INVALID_WORK_DIMENSION:
          return "CL_INVALID_WORK_DIMENSION";
     case CL_INVALID_WORK_GROUP_SIZE:
          return "CL_INVALID_WORK_GROUP_SIZE";
     case CL_INVALID_WORK_ITEM_SIZE:
          return "CL_INVALID_WORK_ITEM_SIZE";
     case CL_INVALID_GLOBAL_OFFSET:
          return "CL_INVALID_GLOBAL_OFFSET";
     case CL_INVALID_EVENT_WAIT_LIST:
          return "CL_INVALID_EVENT_WAIT_LIST";
     case CL_INVALID_EVENT:
          return "CL_INVALID_EVENT";
     case CL_INVALID_OPERATION:
          return "CL_INVALID_OPERATION";
     case CL_INVALID_GL_OBJECT:
          return "CL_INVALID_GL_OBJECT";
     case CL_INVALID_BUFFER_SIZE:
          return "CL_INVALID_BUFFER_SIZE";
     case CL_INVALID_MIP_LEVEL:
          return "CL_INVALID_MIP_LEVEL";
     case CL_INVALID_GLOBAL_WORK_SIZE:
          return "CL_INVALID_GLOBAL_WORK_SIZE";
     case CL_INVALID_PROPERTY:
          return "CL_INVALID_PROPERTY";
     case CL_INVALID_IMAGE_DESCRIPTOR:
          return "CL_INVALID_IMAGE_DESCRIPTOR";
     case CL_INVALID_COMPILER_OPTIONS:
          return "CL_INVALID_COMPILER_OPTIONS";
     case CL_INVALID_LINKER_OPTIONS:
          return "CL_INVALID_LINKER_OPTIONS";
     case CL_INVALID_DEVICE_PARTITION_COUNT:
          return "CL_INVALID_DEVICE_PARTITION_COUNT";
     default:
          return "Unknown OpenCL error";
     }
}

void print_error(cl_int res){
     if (res != CL_SUCCESS){
          printf("error %s\n", get_error(res));
     }
}

void print_build_log_failure(cl_int res, cl_device_id device, cl_program program){
     if (res == CL_BUILD_PROGRAM_FAILURE){
          size_t log_size;
          clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
          char *log = (char *)malloc(log_size);
          clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
          printf("%s\n", log);
     }
}

cl_uint get_num_compute_units(cl_device_id device){
     cl_uint num_compute_units;
     clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &num_compute_units, NULL);
     return num_compute_units;
}

size_t get_work_group_size(cl_kernel kernel, cl_device_id device){
     size_t work_group_size;
     clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
     return work_group_size;
}

cl_uint get_num_dimensions(cl_device_id device){
     cl_uint num_dimensions;
     clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &num_dimensions, NULL);
     return num_dimensions;
}

size_t* get_max_num_work_item(cl_device_id device, cl_uint num_dimensions){
     size_t n = sizeof(size_t) * num_dimensions;
     size_t* max_num_work_item = (size_t*) malloc(n);
     clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, n, max_num_work_item, NULL);
     return max_num_work_item;
}

size_t* get_local_work_size(cl_kernel kernel, cl_device_id device){
     size_t* max_num_work_item = get_max_num_work_item(device, get_num_dimensions(device));
     size_t work_group_size = get_work_group_size(kernel, device);
     size_t* local_work_size = (size_t*) malloc(sizeof(size_t) * 2);
     local_work_size[0] = MIN(floor(sqrt(work_group_size)), max_num_work_item[0]);
     local_work_size[1] = MIN(floor(sqrt(work_group_size)), max_num_work_item[1]);
     return local_work_size;
}

size_t* get_global_work_size(size_t* local_work_size, int width, int height){
     size_t* global_work_size = (size_t*) malloc(sizeof(size_t) * 2);
     global_work_size[0] = ceil(width/(cl_double)local_work_size[0]) * local_work_size[0];
     global_work_size[1] = ceil(height/(cl_double)local_work_size[1]) * local_work_size[1];
     return global_work_size;
}

double* run_kernel(
     char* kernel_filename,
     char* function,
     void* host_ptr1,
     void* host_ptr2,
     void* host_ptr3,
     size_t size1,
     size_t size2,
     size_t size3,
     unsigned int width,
     unsigned int height){
     
     cl_context context;
     cl_command_queue queue;
     cl_program program;
     cl_kernel kernel;
     cl_int errcode_ret;
     cl_mem device_par1_buffer;
     cl_mem device_par2_buffer;
     cl_mem device_par3_buffer;
     cl_mem device_output_buffer;
     size_t output_size;
     double *output_buffer;
     cl_uint num_compute_units;
     char *kernel_source;
     cl_int res;
     struct timeb start, end;
     cl_platform_id gpu_platform;
     cl_device_id device;

     ftime(&start);

     get_GPU_platform(&gpu_platform);
     get_GPU_device(gpu_platform, &device);

     output_size = width * height * sizeof(double);
     output_buffer = (double *)malloc(output_size);
     num_compute_units = get_num_compute_units(device);
     kernel_source = read_file(kernel_filename);     

     context = clCreateContext(0, 1, &device, NULL, NULL, &errcode_ret);
     queue = clCreateCommandQueue(context, device, 0, &errcode_ret);
     program = clCreateProgramWithSource(context, 1, &kernel_source, NULL, &errcode_ret);
     res = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
     print_build_log_failure(res, device, program);
     kernel = clCreateKernel(program, function, &errcode_ret);     
     
     device_par1_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size1, host_ptr1, NULL);
     //device_OpenSimplexEnv_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(OpenSimplexEnv), NULL, NULL);
     //clEnqueueWriteBuffer(queue, device_OpenSimplexEnv_buffer, CL_TRUE, 0, sizeof(OpenSimplexEnv), &ose, 0, NULL, NULL);
     device_par2_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size2, host_ptr2, NULL);
     //device_OpenSimplexGradients_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(OpenSimplexGradients), NULL, NULL);
     //clEnqueueWriteBuffer(queue, device_OpenSimplexGradients_buffer, CL_TRUE, 0, sizeof(OpenSimplexGradients), &osg, 0, NULL, NULL);
     device_par3_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size3, host_ptr3, NULL);
     device_output_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, output_size, NULL, NULL);
     //device_output_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, output_size, output_buffer, NULL);

     clSetKernelArg(kernel, 0, sizeof(cl_mem), &device_par1_buffer);
     clSetKernelArg(kernel, 1, sizeof(cl_mem), &device_par2_buffer);
     clSetKernelArg(kernel, 2, sizeof(cl_mem), &device_par3_buffer);
     clSetKernelArg(kernel, 3, sizeof(unsigned int), &width);
     clSetKernelArg(kernel, 4, sizeof(unsigned int), &height);
     clSetKernelArg(kernel, 5, sizeof(cl_mem), &device_output_buffer);

     size_t* local_work_size = get_local_work_size(kernel, device);
     size_t* global_work_size = get_global_work_size(local_work_size, width, height);

     
     res = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL);
     print_error(res);
     clFinish(queue);
     
     clEnqueueReadBuffer(queue, device_output_buffer, CL_TRUE, 0, output_size, output_buffer, 0, NULL, NULL);
     //output_buffer = (double *) clEnqueueMapBuffer(queue, device_output_buffer, CL_TRUE, CL_MAP_READ, 0, output_size, 0, NULL, NULL, &res);

     /*if (output_buffer == NULL){
          print_error(res);
     }*/

     clReleaseMemObject(device_par1_buffer);
     clReleaseMemObject(device_par2_buffer);
     clReleaseMemObject(device_par3_buffer);
     clReleaseMemObject(device_output_buffer);
     clReleaseProgram(program);
     clReleaseKernel(kernel);
     clReleaseCommandQueue(queue);
     clReleaseContext(context);

     ftime(&end);
     printf("time %s: %fs\n", function, get_time_s(start, end));

     return output_buffer;
}

double* generate_noise2(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise2",
          osg->perm,
          osg->permGrad2,
          ose->LOOKUP_2D,
          sizeof(osg->perm),
          sizeof(osg->permGrad2),
          sizeof(ose->LOOKUP_2D),
          WIDTH,
          HEIGHT);
}

double* generate_noise2_XBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise2_XBeforeY",
          osg->perm,
          osg->permGrad2,
          ose->LOOKUP_2D,
          sizeof(osg->perm),
          sizeof(osg->permGrad2),
          sizeof(ose->LOOKUP_2D),
          WIDTH,
          HEIGHT);
}

double* generate_noise3_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise3_Classic",
          osg->perm,
          osg->permGrad3,
          ose->LOOKUP_3D,
          sizeof(osg->perm),
          sizeof(osg->permGrad3),
          sizeof(ose->LOOKUP_3D),
          WIDTH,
          HEIGHT);
}

double* generate_noise3_XYBeforeZ(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise3_XYBeforeZ",
          osg->perm,
          osg->permGrad3,
          ose->LOOKUP_3D,
          sizeof(osg->perm),
          sizeof(osg->permGrad3),
          sizeof(ose->LOOKUP_3D),
          WIDTH,
          HEIGHT);
}

double* generate_noise3_XZBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise3_XZBeforeY",
          osg->perm,
          osg->permGrad3,
          ose->LOOKUP_3D,
          sizeof(osg->perm),
          sizeof(osg->permGrad3),
          sizeof(ose->LOOKUP_3D),
          WIDTH,
          HEIGHT);
}

double* generate_noise4_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise4_Classic",
          osg->perm,
          osg->permGrad4,
          ose->VERTICES_4D,
          sizeof(osg->perm),
          sizeof(osg->permGrad4),
          sizeof(ose->VERTICES_4D),
          WIDTH,
          HEIGHT);
}

double* generate_noise4_XYBeforeZW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise4_XYBeforeZW",
          osg->perm,
          osg->permGrad4,
          ose->VERTICES_4D,
          sizeof(osg->perm),
          sizeof(osg->permGrad4),
          sizeof(ose->VERTICES_4D),
          WIDTH,
          HEIGHT);
}

double* generate_noise4_XZBeforeYW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise4_XZBeforeYW",
          osg->perm,
          osg->permGrad4,
          ose->VERTICES_4D,
          sizeof(osg->perm),
          sizeof(osg->permGrad4),
          sizeof(ose->VERTICES_4D),
          WIDTH,
          HEIGHT);
}

double* generate_noise4_XYZBeforeW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
     return run_kernel(
          "OpenSimplex2F.cl",
          "noise4_XYZBeforeW",
          osg->perm,
          osg->permGrad4,
          ose->VERTICES_4D,
          sizeof(osg->perm),
          sizeof(osg->permGrad4),
          sizeof(ose->VERTICES_4D),
          WIDTH,
          HEIGHT);
}

int main(){
     OpenSimplexEnv ose = initOpenSimplex();
     OpenSimplexGradients osg = newOpenSimplexGradients(&ose, 1234);

     save_bitmap("img/noise2.bmp", WIDTH, HEIGHT, generate_noise2(&ose, &osg));
     save_bitmap("img/noise2_XBeforeY.bmp", WIDTH, HEIGHT, generate_noise2_XBeforeY(&ose, &osg));
     save_bitmap("img/noise3_Classic.bmp", WIDTH, HEIGHT, generate_noise3_Classic(&ose, &osg));
     save_bitmap("img/noise3_XYBeforeZ.bmp", WIDTH, HEIGHT, generate_noise3_XYBeforeZ(&ose, &osg));
     save_bitmap("img/noise3_XZBeforeY.bmp", WIDTH, HEIGHT, generate_noise3_XZBeforeY(&ose, &osg));
     save_bitmap("img/noise4_Classic.bmp", WIDTH, HEIGHT, generate_noise4_Classic(&ose, &osg));
     save_bitmap("img/noise4_XYBeforeZW.bmp", WIDTH, HEIGHT, generate_noise4_XYBeforeZW(&ose, &osg));
     save_bitmap("img/noise4_XZBeforeYW.bmp", WIDTH, HEIGHT, generate_noise4_XZBeforeYW(&ose, &osg));
     save_bitmap("img/noise4_XYZBeforeW.bmp", WIDTH, HEIGHT, generate_noise4_XYZBeforeW(&ose, &osg));
     return 0;
}