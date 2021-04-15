#include <stdio.h>
#include "../OpenSimplex2F.h"



int main(){
    OpenSimplexEnv *ose = initOpenSimplex();
    OpenSimplexGradients *osg = newOpenSimplexGradients(ose, 1234);
    for (int x = 0; x < 10; x++){
        for (int y = 0; y < 10; y++){
            double val = noise2(ose, osg, x, y);
            printf("val %f\n", val);
        }
    }

    for (int x = 0; x < 10; x++){
        for (int y = 0; y < 10; y++){
            for (int z = 0; z < 10; z++){
                double val = noise3_Classic(ose, osg, x, y, z);
                printf("val %f\n", val);
            }
        }
    }

    return 0;
}