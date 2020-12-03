// C++ headers
#include "headcpp.h"

// C headers
#include "stdlib.h"
#include <stdio.h>
#include <string.h>
#include "assert.h"
#include "math.h"

// Lorene headers
#include "nbr_spx.h"
#include "tensor.h"
#include "metric.h"
#include "graphique.h"
#include "utilitaires.h"

/*
TODO:
    - Use get_bvect_cart() since this is read as cartesian!
      Then use change_triad(map.get_bvect_spher())

*/

using namespace Lorene ;



int main(int argc, char **argv) {
    // Reads parameter file.
    // Could be done in a function, but that is a hasle...



    FILE* fp = fopen("lorene_parameters.txt", "r");
    int line = 0;
    char buf[255];
    
    int nz;

    if (fgets(buf, sizeof(buf), fp) != NULL){
        nz = atoi(buf);
        cout << nz << endl;
    }


    int nr_array[nz];
    int nt_array[nz];
    int np_array[nz];
    int type_r[nz];
    double r_limits[nz+2];

    double temp_holder[nz*4+2];
    int temp_countr = 0;

    while( fgets(buf, sizeof(buf), fp) != NULL){
        char *tolken = strtok(buf, "\n"); 
        //cout << tolken << endl;
        char *token = strtok(tolken, " ");
        while(token != NULL){
            temp_holder[temp_countr++] = atof(token);
            token = strtok(NULL, " ");
        }
    }

    for(int i = 0; i<nz; i++){
        nr_array[i] = temp_holder[i];
    }
    for(int i = 0; i<nz; i++){
        nt_array[i] = temp_holder[nz+i];
    }
    for(int i = 0; i<nz; i++){
        np_array[i] = temp_holder[2*nz+i];
    }

    for(int i = 0; i<nz+2; i++){
        r_limits[i] = temp_holder[3*nz+i];
    }

    r_limits[nz+1] = __infinity;

    type_r[0] = RARE;
    type_r[nz-1] = UNSURR;
    for(int i = 1; i < nz-1; i++){
        type_r[i] = FIN;
    }

 
    








}