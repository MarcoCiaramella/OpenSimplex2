#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../OpenSimplex2F/OpenSimplex2F.h"



#define SIZE 512
#define OFF_X 2048
#define OFF_Y 2048
#define OFF_Z 2048
#define FREQ1 1.0
#define FREQ2 2.0
#define FREQ4 4.0
#define FREQ8 8.0
#define FREQ16 16.0
#define FREQ32 32.0
#define OCT1 0.9
#define OCT2 0.6
#define OCT4 0.5
#define OCT8 0.3
#define OCT16 0.2
#define OCT32 0.1
#define EXP 3.0


typedef struct {
    float* vertices;
    float* normals;
    int* indices;
    int num_vertices;
    int num_triangles;
} mesh;




mesh add_noise(OpenSimplexEnv *ose, OpenSimplexGradients *osg, mesh m, double (*noise_fun)(OpenSimplexEnv*, OpenSimplexGradients*, double, double, double)){
    for (int i = 0; i < m.num_vertices; i++){
        int j = i*3;

        float vx = m.vertices[j];
        float vy = m.vertices[j+1];
        float vz = m.vertices[j+2];

        float noise1 = noise_fun(ose, osg, (vx + OFF_X) * FREQ1, (vy + OFF_Y) * FREQ1, (vz + OFF_Z) * FREQ1);
        float noise2 = noise_fun(ose, osg, (vx + OFF_X) * FREQ2, (vy + OFF_Y) * FREQ2, (vz + OFF_Z) * FREQ2);
        float noise4 = noise_fun(ose, osg, (vx + OFF_X) * FREQ4, (vy + OFF_Y) * FREQ4, (vz + OFF_Z) * FREQ4);
        float noise8 = noise_fun(ose, osg, (vx + OFF_X) * FREQ8, (vy + OFF_Y) * FREQ8, (vz + OFF_Z) * FREQ8);
        float noise16 = noise_fun(ose, osg, (vx + OFF_X) * FREQ16, (vy + OFF_Y) * FREQ16, (vz + OFF_Z) * FREQ16);
        float noise32 = noise_fun(ose, osg, (vx + OFF_X) * FREQ32, (vy + OFF_Y) * FREQ32, (vz + OFF_Z) * FREQ32);

        float noise = (OCT1*noise1 + OCT2*noise2 + OCT4*noise4 + OCT8*noise8 + OCT16*noise16 + OCT32*noise32 + 1.0) / 2.0;
        noise /= OCT1 + OCT2 + OCT4 + OCT8 + OCT16 + OCT32;
        noise = pow(noise, EXP);

        float nx = m.normals[j];
        float ny = m.normals[j+1];
        float nz = m.normals[j+2];
        m.vertices[j] = vx + nx*noise;
        m.vertices[j+1] = vy + ny*noise;
        m.vertices[j+2] = vz + nz*noise;
    }
    return m;
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
    export_mesh_ply("planet_noise3_Classic.ply", add_noise(ose, osg, new_sphere(2, SIZE, SIZE), noise3_Classic));
    export_mesh_ply("planet_noise3_XYBeforeZ.ply", add_noise(ose, osg, new_sphere(2, SIZE, SIZE), noise3_XYBeforeZ));
    export_mesh_ply("planet_noise3_XZBeforeY.ply", add_noise(ose, osg, new_sphere(2, SIZE, SIZE), noise3_XZBeforeY));
    return 0;
}