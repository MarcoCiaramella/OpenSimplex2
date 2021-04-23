#include <stdlib.h>
#include <stdio.h>
#include "../OpenSimplex2F.h"
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

double **generate_noise3(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise3_Classic(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ);
        }
    }
    return noise;
}

double **generate_noise4(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise4_Classic(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ, 0.0);
        }
    }
    return noise;
}

void generate_testing_images(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    save_bitmap("test/img/noise2.bmp", WIDTH, HEIGHT, generate_noise2(ose, osg));
    save_bitmap("test/img/noise3.bmp", WIDTH, HEIGHT, generate_noise3(ose, osg));
    save_bitmap("test/img/noise4.bmp", WIDTH, HEIGHT, generate_noise4(ose, osg));
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

void test_noise3(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise3.bmp", save_bitmap("test/img/noise3_tmp.bmp", WIDTH, HEIGHT, generate_noise3(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise3 passed: %fs\n", s);
    }
    else {
        printf("test noise3 not passed: %fs\n", s);
    }
    remove("test/img/noise3_tmp.bmp");
}

void test_noise4(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    struct timeb start, end;
    ftime(&start);
    int res = bitmap_compare("test/img/noise4.bmp", save_bitmap("test/img/noise4_tmp.bmp", WIDTH, HEIGHT, generate_noise4(ose, osg)));
    ftime(&end);
    float s = get_time_s(start, end);
    if (res){
        printf("test noise4 passed: %fs\n", s);
    }
    else {
        printf("test noise4 not passed: %fs\n", s);
    }
    remove("test/img/noise4_tmp.bmp");
}


int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    test_noise2(ose, osg);
    test_noise3(ose, osg);
    test_noise4(ose, osg);
    return 0;
}