# OpenSimplex2
CPU implementation in C of [OpenSimplex 2](https://github.com/KdotJPG/OpenSimplex2)

## How to use
```c
#include <stdlib.h>
#include <stdio.h>
// use #include "../OpenSimplex2S.h" for OpenSimplex2S
#include "../OpenSimplex2F.h"


#define BYTES_PER_PIXEL 3
#define FILE_HEADER_SIZE 14
#define INFO_HEADER_SIZE 40



unsigned char *createBitmapInfoHeader(int height, int width){
    static unsigned char infoHeader[] = {
        0,0,0,0, /// header size
        0,0,0,0, /// image width
        0,0,0,0, /// image height
        0,0,     /// number of color planes
        0,0,     /// bits per pixel
        0,0,0,0, /// compression
        0,0,0,0, /// image size
        0,0,0,0, /// horizontal resolution
        0,0,0,0, /// vertical resolution
        0,0,0,0, /// colors in color table
        0,0,0,0, /// important color count
    };

    infoHeader[ 0] = (unsigned char)(INFO_HEADER_SIZE);
    infoHeader[ 4] = (unsigned char)(width      );
    infoHeader[ 5] = (unsigned char)(width >>  8);
    infoHeader[ 6] = (unsigned char)(width >> 16);
    infoHeader[ 7] = (unsigned char)(width >> 24);
    infoHeader[ 8] = (unsigned char)(height      );
    infoHeader[ 9] = (unsigned char)(height >>  8);
    infoHeader[10] = (unsigned char)(height >> 16);
    infoHeader[11] = (unsigned char)(height >> 24);
    infoHeader[12] = (unsigned char)(1);
    infoHeader[14] = (unsigned char)(BYTES_PER_PIXEL*8);

    return infoHeader;
}

unsigned char *createBitmapFileHeader(int height, int stride){
    int fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

    static unsigned char fileHeader[] = {
        0,0,     /// signature
        0,0,0,0, /// image file size in bytes
        0,0,0,0, /// reserved
        0,0,0,0, /// start of pixel array
    };

    fileHeader[ 0] = (unsigned char)('B');
    fileHeader[ 1] = (unsigned char)('M');
    fileHeader[ 2] = (unsigned char)(fileSize      );
    fileHeader[ 3] = (unsigned char)(fileSize >>  8);
    fileHeader[ 4] = (unsigned char)(fileSize >> 16);
    fileHeader[ 5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = (unsigned char)(FILE_HEADER_SIZE + INFO_HEADER_SIZE);

    return fileHeader;
}

void generateBitmapImage(unsigned char *image, int height, int width, char *imageFileName){
    int widthInBytes = width * BYTES_PER_PIXEL;

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;

    int stride = (widthInBytes) + paddingSize;

    FILE* imageFile = fopen(imageFileName, "wb");

    unsigned char* fileHeader = createBitmapFileHeader(height, stride);
    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);

    unsigned char* infoHeader = createBitmapInfoHeader(height, width);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

    int i;
    for (i = 0; i < height; i++) {
        fwrite(image + (i*widthInBytes), BYTES_PER_PIXEL, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
}

void save_bitmap(char *filename, int width, int height, double **vals){
    unsigned char *image = (unsigned char *) malloc(sizeof(unsigned char) * height * width * BYTES_PER_PIXEL);
    int index = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            double val = vals[i][j];
            val = (val + 1.0) / 2.0;
            unsigned char gray = (unsigned char) (val * 255);
            image[index++] = gray;
            image[index++] = gray;
            image[index++] = gray;
        }
    }

    generateBitmapImage(image, height, width, filename);
    free(image);
}



#define WIDTH 512
#define HEIGHT 512
#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD



int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = noise2(ose, osg, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ);
        }
    }
    save_bitmap("noise2.bmp", WIDTH, HEIGHT, noise);
    free(noise);
    return 0;
}
```

## Performance
Image 4096x4096 on CPU Intel Core i5-4460 3.20GHz
### OpenSimplex2F
* noise2 1.793000s
* noise2_XBeforeY 1.733000s
* noise3_Classic 2.108000s
* noise3_XYBeforeZ 2.111000s
* noise3_XZBeforeY 2.149000s
* noise4_Classic 3.126000s
* noise4_XYBeforeZW 3.423000s
* noise4_XZBeforeYW 3.234000s
* noise4_XYZBeforeW 3.155000s
### OpenSimplex2S
* noise2 2.162000s
* noise2_XBeforeY 2.092000s
* noise3_Classic 2.854000s
* noise3_XYBeforeZ 2.861000s
* noise3_XZBeforeY 2.851000s
* noise4_Classic 5.153000s
* noise4_XYBeforeZW 5.397000s
* noise4_XZBeforeYW 5.312000s
* noise4_XYZBeforeW 4.971000s