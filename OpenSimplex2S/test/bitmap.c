#include <stdlib.h>
#include <stdio.h>



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

void generateBitmapImage(unsigned char *image, int height, int width, FILE *imageFile){
    int widthInBytes = width * BYTES_PER_PIXEL;

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;

    int stride = (widthInBytes) + paddingSize;

    unsigned char* fileHeader = createBitmapFileHeader(height, stride);
    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);

    unsigned char* infoHeader = createBitmapInfoHeader(height, width);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

    int i;
    for (i = 0; i < height; i++) {
        fwrite(image + (i*widthInBytes), BYTES_PER_PIXEL, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }
}

unsigned char *create_bitmap(int width, int height, double **vals){
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
    return image;
}

char *save_bitmap(char *filename, int width, int height, double **vals){
    unsigned char *image = create_bitmap(width, height, vals);
    FILE *f = fopen(filename, "wb");
    generateBitmapImage(image, height, width, f);
    fclose(f);
    free(image);
    return filename;
}

unsigned char *load_bitmap(char *filename, long *num_bytes){
    FILE *f = fopen(filename,"rb");
    unsigned char info[54];

    // read the 54-byte header
    fread(info, sizeof(unsigned char), 54, f); 

    // extract image height and width from header
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];

    // allocate 3 bytes per pixel
    int size = BYTES_PER_PIXEL * width * height;
    unsigned char *data = (unsigned char *) malloc(sizeof(unsigned char) * size);

    // read the rest of the data at once
    fread(data, sizeof(unsigned char), size, f); 
    fclose(f);

    for (int i = 0; i < size; i += BYTES_PER_PIXEL){
        // flip the order of every 3 bytes
        unsigned char tmp = data[i];
        data[i] = data[i+2];
        data[i+2] = tmp;
    }

    *num_bytes = sizeof(unsigned char) * size;

    return data;
}

int bitmap_compare(char *filename1, char *filename2){
    long num_bytes;
    unsigned char *data1 = load_bitmap(filename1, &num_bytes);
    unsigned char *data2 = load_bitmap(filename2, &num_bytes);
    return memcmp(data1, data2, num_bytes) == 0;
}