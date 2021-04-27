#include <stdio.h>
#include <CL/cl.h>



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

void get_GPU_device(cl_platform_id platform, cl_device_id* gpu_device){
     if (get_num_GPU_devices(platform) > 0){
          cl_device_id devices[1];
          exit_on_error(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, devices, NULL));
          *gpu_device = devices[0];
          return;
     }
     printf("No GPU device found.\n");
     exit(-1);
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
          if (strcmp(platform_name, "NVIDIA CUDA") == 0){
               *gpu_platform = platforms[i];
               return;
          }
    }
    printf("No GPU platform found.\n");
    exit(-1);
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

void print_build_log_failure(cl_int res, cl_device_id gpu_device, cl_program program){
     if (res == CL_BUILD_PROGRAM_FAILURE){
          size_t log_size;
          clGetProgramBuildInfo(program, gpu_device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
          char *log = (char *)malloc(log_size);
          clGetProgramBuildInfo(program, gpu_device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
          printf("%s\n", log);
     }
}

void run_kernel(cl_device_id gpu_device, char *kernel_filename, int width, int height){
     cl_context context;
     cl_command_queue queue;
     cl_program program;
     cl_kernel kernel;
     cl_int errcode_ret;
     cl_mem device_output_buffer;
     double *output_buffer;
     unsigned int size = width * height;
     size_t buffer_size_in_bytes = size * sizeof(double);
     size_t num_work_groups[] = {width, height};
     size_t work_group_size[] = {1, 1};
     output_buffer = (double *)malloc(buffer_size_in_bytes);
     char *kernel_source = read_file(kernel_filename);

     context = clCreateContext(0, 1, &gpu_device, NULL, NULL, &errcode_ret);
     queue = clCreateCommandQueue(context, gpu_device, 0, &errcode_ret);
     program = clCreateProgramWithSource(context, 1, &kernel_source, NULL, &errcode_ret);
     cl_int res = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
     print_build_log_failure(res, gpu_device, program);
     kernel = clCreateKernel(program, "main", &errcode_ret);
     device_output_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, buffer_size_in_bytes, NULL, NULL);

     clSetKernelArg(kernel, 0, sizeof(unsigned int), &size);
     clSetKernelArg(kernel, 1, sizeof(cl_mem), &device_output_buffer);

     clEnqueueNDRangeKernel(queue, kernel, 2, NULL, num_work_groups, work_group_size, 0, NULL, NULL);

     clFinish(queue);

     clEnqueueReadBuffer(queue, device_output_buffer, CL_TRUE, 0, buffer_size_in_bytes, output_buffer, 0, NULL, NULL);

     for (int i = 0; i < size; i++){
          printf("val %f\n", output_buffer[i]);
     }

     clReleaseMemObject(device_output_buffer);
     clReleaseProgram(program);
     clReleaseKernel(kernel);
     clReleaseCommandQueue(queue);
     clReleaseContext(context);
     free(output_buffer);
}

int main(){
     cl_platform_id gpu_platform;
     cl_device_id gpu_device;
     get_GPU_platform(&gpu_platform);
     get_GPU_device(gpu_platform, &gpu_device);
     run_kernel(gpu_device, "OpenSimplex2F.cl", 10, 10);
     return 0;
}