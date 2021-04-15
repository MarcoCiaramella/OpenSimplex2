#include <stdlib.h>
#include <stdio.h>
#include "../OpenSimplex2F.h"



#define WIDTH 100
#define HEIGHT 100



unsigned char *to_rgb_arr(int width, int height, float **vals){
    int size = width * height * 4;
    unsigned char *pixels = (unsigned char *) malloc(size);
    for (int x = 0; x < width; x++){
        for (int y = 0; y < height; y++){
            float val = vals[x][y];
            val = (val + 1.0) / 2.0;
            unsigned char gray = (unsigned char) (val * 255);
            int p = (y * width + x) * 4;
            pixels[p + 0] = gray; //blue
            pixels[p + 1] = gray;//green
            pixels[p + 2] = gray;//red
        }
    }
    return pixels;
}

void save_bitmap(char *filename, int width, int height, float **vals){
    int size = width * height * 4; //for 32-bit bitmap only

    char header[54] = { 0 };
    strcpy(header, "BM");
    memset(&header[2],  (int)(54 + size), 1);
    memset(&header[10], (int)54, 1);//always 54
    memset(&header[14], (int)40, 1);//always 40
    memset(&header[18], (int)width, 1);
    memset(&header[22], (int)height, 1);
    memset(&header[26], (short)1, 1);
    memset(&header[28], (short)32, 1);//32bit
    memset(&header[34], (int)size, 1);//pixel size

    unsigned char *pixels = to_rgb_arr(width, height, vals);

    FILE *fout = fopen(filename, "wb");
    fwrite(header, 1, 54, fout);
    fwrite(pixels, 1, size, fout);
    free(pixels);
    fclose(fout);
}

int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    float **noise = (float **) malloc(sizeof(float *) * WIDTH);
    for (int x = 0; x < WIDTH; x++){
        noise[x] = (float *) malloc(sizeof(float) * HEIGHT);
        for (int y = 0; y < HEIGHT; y++){
            noise[x][y] = noise2(ose, osg, x, y);
        }
    }
    save_bitmap("noise2.bmp", WIDTH, HEIGHT, noise);

    for (int x = 0; x < 10; x++){
        for (int y = 0; y < 10; y++){
            for (int z = 0; z < 10; z++){
                double val = noise3_Classic(ose, osg, x, y, z);
            }
        }
    }

    return 0;
}