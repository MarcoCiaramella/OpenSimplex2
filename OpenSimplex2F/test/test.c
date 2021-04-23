#include <stdlib.h>
#include <stdio.h>
#include "../OpenSimplex2F.h"
#include "bitmap.h"



#define WIDTH 512
#define HEIGHT 512
#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD




void generate_testing_images(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise2(ose, osg, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ);
        }
    }
    save_bitmap("test/img/noise2.bmp", WIDTH, HEIGHT, noise);
    free(noise);

    noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise3_Classic(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ);
        }
    }
    save_bitmap("test/img/noise3.bmp", WIDTH, HEIGHT, noise);
    free(noise);
}

void test_noise2(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise2(ose, osg, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ);
        }
    }
    if (bitmap_compare("test/img/noise2.bmp", save_bitmap("test/img/noise2_tmp.bmp", WIDTH, HEIGHT, noise))){
        printf("ok\n");
    }
    else {
        printf("not ok\n");
    }
    remove("test/img/noise2_tmp.bmp");
    free(noise);
}

void test_noise3(OpenSimplexEnv *ose, OpenSimplexGradients *osg){
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise3_Classic(ose, osg, (x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ);
        }
    }
    if (bitmap_compare("test/img/noise3.bmp", save_bitmap("test/img/noise3_tmp.bmp", WIDTH, HEIGHT, noise))){
        printf("ok\n");
    }
    else {
        printf("not ok\n");
    }
    remove("test/img/noise3_tmp.bmp");
    free(noise);
}



int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    test_noise2(ose, osg);
    test_noise3(ose, osg);
    return 0;
}