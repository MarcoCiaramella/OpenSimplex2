#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../OpenSimplex2F.h"



#define WIDTH 4096
#define HEIGHT 4096
#define PERIOD 1.0
#define OFF_X 2048
#define OFF_Y 2048
#define OFF_Z 2048
#define FREQ 1.0 / PERIOD


typedef struct {
    float* vertices;
    float* normals;
    int* indices;
    int num_vertices;
    int num_triangles;
} mesh;




void generate_noise3_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg, mesh m){
    for (int i = 0; i < m.num_vertices; i++){
        int j = i*3;
        float vx = m.vertices[j];
        float vy = m.vertices[j+1];
        float vz = m.vertices[j+2];
        float noise = noise3_Classic(ose, osg, (vx + OFF_X) * FREQ, (vy + OFF_Y) * FREQ, (vz + OFF_Z) * FREQ);
        noise = (noise+1.0)/2.0;
        //float noise = noise3_Classic(ose, osg, (vx + OFF_X) * FREQ, 0.0, (vy + OFF_Y) * FREQ);
        float nx = m.normals[j];
        float ny = m.normals[j+1];
        float nz = m.normals[j+2];
        m.vertices[j] = vx + nx*noise;
        m.vertices[j+1] = vy + ny*noise;
        m.vertices[j+2] = vz + nz*noise;
    }
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

// sectorCount;                        // longitude, # of slices
// stackCount;                         // latitude, # of stacks
mesh new_sphere(float radius, int sectorCount, int stackCount){
    const float PI = acos(-1);

    mesh sphere;
    sphere.vertices = (float*) malloc(sizeof(float) * (sectorCount+1) * (stackCount+1) * 3);
    sphere.normals = (float*) malloc(sizeof(float) * (sectorCount+1) * (stackCount+1) * 3);
    sphere.indices = (int*) malloc(sizeof(int) * (sectorCount + sectorCount + (stackCount-2)*sectorCount*2) * 3);

    float x, y, z, xy;                              // vertex position
    float nx, ny, nz, lengthInv = 1.0f / radius;    // normal

    float sectorStep = 2 * PI / sectorCount;
    float stackStep = PI / stackCount;
    float sectorAngle, stackAngle;

    int iv = 0;
    int in = 0;
    int ii = 0;

    sphere.num_vertices = 0;
    sphere.num_triangles = 0;

    for(int i = 0; i <= stackCount; ++i){
        stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = radius * cosf(stackAngle);             // r * cos(u)
        z = radius * sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for(int j = 0; j <= sectorCount; ++j){
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position
            x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
            sphere.vertices[iv++] = x;
            sphere.vertices[iv++] = y;
            sphere.vertices[iv++] = z;

            // normalized vertex normal
            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;
            sphere.normals[in++] = nx;
            sphere.normals[in++] = ny;
            sphere.normals[in++] = nz;

            sphere.num_vertices++;
        }
    }

    // indices
    //  k1--k1+1
    //  |  / |
    //  | /  |
    //  k2--k2+1
    unsigned int k1, k2;
    for(int i = 0; i < stackCount; ++i){
        k1 = i * (sectorCount + 1);     // beginning of current stack
        k2 = k1 + sectorCount + 1;      // beginning of next stack

        for(int j = 0; j < sectorCount; ++j, ++k1, ++k2){
            // 2 triangles per sector excluding 1st and last stacks
            if(i != 0){
                // k1---k2---k1+1
                sphere.indices[ii++] = k1;
                sphere.indices[ii++] = k2;
                sphere.indices[ii++] = k1+1;

                sphere.num_triangles++;
            }

            if(i != (stackCount-1)){
                // k1+1---k2---k2+1
                sphere.indices[ii++] = k1+1;
                sphere.indices[ii++] = k2;
                sphere.indices[ii++] = k2+1;

                sphere.num_triangles++;
            }
        }
    }

    return sphere;
}

void export_mesh_ply(const char* ply_filename, const mesh m){

    const char* header = "ply\n"
    "format ascii 1.0\n"
    "element vertex %d\n"
    "property float x\n"
    "property float y\n"
    "property float z\n"
    "property float nx\n"
    "property float ny\n"
    "property float nz\n"
    "element face %d\n"
    "property list uchar uint vertex_indices\n"
    "end_header\n";
    const char* vertex = "%f %f %f %f %f %f\n";
    const char* face = "3 %d %d %d\n";


    FILE* f = fopen(ply_filename, "w");
    char buffer[255];
    int n = sprintf(buffer, header, m.num_vertices, m.num_triangles);
    fwrite(buffer, 1, n, f);
    for (int i = 0; i < m.num_vertices; i++){
        int j = i*3;
        n = sprintf(
            buffer,
            vertex,
            m.vertices[j],
            m.vertices[j+1],
            m.vertices[j+2],
            m.normals[j],
            m.normals[j+1],
            m.normals[j+2]);
        fwrite(buffer, 1, n, f);
    }
    for (int i = 0; i < m.num_triangles; i++){
        int j = i*3;
        n = sprintf(
            buffer,
            face,
            m.indices[j],
            m.indices[j+1],
            m.indices[j+2]);
        fwrite(buffer, 1, n, f);
    }
    fclose(f);
}

int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    mesh sphere = new_sphere(2, 512, 512);
    //export_mesh_ply("sphere.ply", sphere);
    generate_noise3_Classic(ose, osg, sphere);
    export_mesh_ply("sphere_noise.ply", sphere);
    return 0;
}