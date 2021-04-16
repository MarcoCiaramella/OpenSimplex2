#include <stdlib.h>
#include <stdio.h>
#include "../OpenSimplex2S.h"
#include "bitmap.h"



int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    float **noise = (float **) malloc(sizeof(float *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (float *) malloc(sizeof(float) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise2(ose, osg, x, y);
        }
    }
    save_bitmap("noise2.bmp", noise);

    for (int x = 0; x < 10; x++){
        for (int y = 0; y < 10; y++){
            for (int z = 0; z < 10; z++){
                double val = noise3_Classic(ose, osg, x, y, z);
            }
        }
    }

    return 0;
}