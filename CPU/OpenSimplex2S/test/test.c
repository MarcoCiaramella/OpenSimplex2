#include <stdlib.h>
#include <stdio.h>
#include "../OpenSimplex2S.h"
#include "bitmap.h"
#include <sys/timeb.h>



#define WIDTH 512
#define HEIGHT 512
#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD




float get_time_s(struct timeb start, struct timeb end){
    unsigned long long ms = (unsigned long long) (1000.0 * (end.time - start.time) + (end.millitm - start.millitm));
    return ms/1000.0;
}

double **generate_noise2(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise2(ose, osg, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ);
        }
    }
    return noise;
}

double **generate_noise2_XBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise2_XBeforeY(ose, osg, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ);
        }
    }
    return noise;
}

double **generate_noise3_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise3_Classic(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ);
        }
    }
    return noise;
}

double **generate_noise3_XYBeforeZ(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise3_XYBeforeZ(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ);
        }
    }
    return noise;
}

double **generate_noise3_XZBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise3_XZBeforeY(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ);
        }
    }
    return noise;
}

double **generate_noise4_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise4_Classic(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ, 0.0);
        }
    }
    return noise;
}

double **generate_noise4_XYBeforeZW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise4_XYBeforeZW(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ, 0.0);
        }
    }
    return noise;
}

double **generate_noise4_XZBeforeYW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise4_XZBeforeYW(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ, 0.0);
        }
    }
    return noise;
}

double **generate_noise4_XYZBeforeW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise4_XYZBeforeW(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ, 0.0);
        }
    }
    return noise;
}

void generate_testing_images(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    save_bitmap("test/img/noise2.bmp", WIDTH, HEIGHT, generate_noise2(ose, osg));
    save_bitmap("test/img/noise2_XBeforeY.bmp", WIDTH, HEIGHT, generate_noise2_XBeforeY(ose, osg));
    save_bitmap("test/img/noise3_Classic.bmp", WIDTH, HEIGHT, generate_noise3_Classic(ose, osg));
    save_bitmap("test/img/noise3_XYBeforeZ.bmp", WIDTH, HEIGHT, generate_noise3_XYBeforeZ(ose, osg));
    save_bitmap("test/img/noise3_XZBeforeY.bmp", WIDTH, HEIGHT, generate_noise3_XZBeforeY(ose, osg));
    save_bitmap("test/img/noise4_Classic.bmp", WIDTH, HEIGHT, generate_noise4_Classic(ose, osg));
    save_bitmap("test/img/noise4_XYBeforeZW.bmp", WIDTH, HEIGHT, generate_noise4_XYBeforeZW(ose, osg));
    save_bitmap("test/img/noise4_XZBeforeYW.bmp", WIDTH, HEIGHT, generate_noise4_XZBeforeYW(ose, osg));
    save_bitmap("test/img/noise4_XYZBeforeW.bmp", WIDTH, HEIGHT, generate_noise4_XYZBeforeW(ose, osg));
}

void test_noise2(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise2.bmp", save_bitmap("test/img/noise2_tmp.bmp", WIDTH, HEIGHT, generate_noise2(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise2 passed: %fs\n", s);
    }
    else {
        printf("test noise2 not passed: %fs\n", s);
    }
    remove("test/img/noise2_tmp.bmp");
}

void test_noise2_XBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise2_XBeforeY.bmp", save_bitmap("test/img/noise2_XBeforeY_tmp.bmp", WIDTH, HEIGHT, generate_noise2_XBeforeY(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise2_XBeforeY passed: %fs\n", s);
    }
    else {
        printf("test noise2_XBeforeY not passed: %fs\n", s);
    }
    remove("test/img/noise2_XBeforeY_tmp.bmp");
}

void test_noise3_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise3_Classic.bmp", save_bitmap("test/img/noise3_Classic_tmp.bmp", WIDTH, HEIGHT, generate_noise3_Classic(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise3_Classic passed: %fs\n", s);
    }
    else {
        printf("test noise3_Classic not passed: %fs\n", s);
    }
    remove("test/img/noise3_Classic_tmp.bmp");
}

void test_noise3_XYBeforeZ(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise3_XYBeforeZ.bmp", save_bitmap("test/img/noise3_XYBeforeZ_tmp.bmp", WIDTH, HEIGHT, generate_noise3_XYBeforeZ(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise3_XYBeforeZ passed: %fs\n", s);
    }
    else {
        printf("test noise3_XYBeforeZ not passed: %fs\n", s);
    }
    remove("test/img/noise3_XYBeforeZ_tmp.bmp");
}

void test_noise3_XZBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise3_XZBeforeY.bmp", save_bitmap("test/img/noise3_XZBeforeY_tmp.bmp", WIDTH, HEIGHT, generate_noise3_XZBeforeY(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise3_XZBeforeY passed: %fs\n", s);
    }
    else {
        printf("test noise3_XZBeforeY not passed: %fs\n", s);
    }
    remove("test/img/noise3_XZBeforeY_tmp.bmp");
}

void test_noise4_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise4_Classic.bmp", save_bitmap("test/img/noise4_Classic_tmp.bmp", WIDTH, HEIGHT, generate_noise4_Classic(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise4_Classic passed: %fs\n", s);
    }
    else {
        printf("test noise4_Classic not passed: %fs\n", s);
    }
    remove("test/img/noise4_Classic_tmp.bmp");
}

void test_noise4_XYBeforeZW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise4_XYBeforeZW.bmp", save_bitmap("test/img/noise4_XYBeforeZW_tmp.bmp", WIDTH, HEIGHT, generate_noise4_XYBeforeZW(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise4_XYBeforeZW passed: %fs\n", s);
    }
    else {
        printf("test noise4_XYBeforeZW not passed: %fs\n", s);
    }
    remove("test/img/noise4_XYBeforeZW_tmp.bmp");
}

void test_noise4_XZBeforeYW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise4_XZBeforeYW.bmp", save_bitmap("test/img/noise4_XZBeforeYW_tmp.bmp", WIDTH, HEIGHT, generate_noise4_XZBeforeYW(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise4_XZBeforeYW passed: %fs\n", s);
    }
    else {
        printf("test noise4_XZBeforeYW not passed: %fs\n", s);
    }
    remove("test/img/noise4_XZBeforeYW_tmp.bmp");
}

void test_noise4_XYZBeforeW(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise4_XYZBeforeW.bmp", save_bitmap("test/img/noise4_XYZBeforeW_tmp.bmp", WIDTH, HEIGHT, generate_noise4_XYZBeforeW(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise4_XYZBeforeW passed: %fs\n", s);
    }
    else {
        printf("test noise4_XYZBeforeW not passed: %fs\n", s);
    }
    remove("test/img/noise4_XYZBeforeW_tmp.bmp");
}


int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    //generate_testing_images(ose, osg);
    test_noise2(ose, osg);
    test_noise2_XBeforeY(ose, osg);
    test_noise3_Classic(ose, osg);
    test_noise3_XYBeforeZ(ose, osg);
    test_noise3_XZBeforeY(ose, osg);
    test_noise4_Classic(ose, osg);
    test_noise4_XYBeforeZW(ose, osg);
    test_noise4_XZBeforeYW(ose, osg);
    test_noise4_XYZBeforeW(ose, osg);
    return 0;
}