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
    
    int dz;

    if (fgets(buf, sizeof(buf), fp) != NULL){
        dz = atoi(buf);
        cout << dz << endl;
    }


    int nr_array[dz];
    int nt_array[dz];
    int np_array[dz];
    int type_r[dz];

    int temp_holder[dz*3];
    int temp_countr = 0;

    while( fgets(buf, sizeof(buf), fp) != NULL){
        char *tolken = strtok(buf, "\n"); 
        cout << tolken << endl;
        char *token = strtok(tolken, " ");
        while(token != NULL){
            temp_holder[temp_countr++] = atoi(token);
            token = strtok(NULL, " ");
        }
    }

    for(int i = 0; i<dz; i++){
        nr_array[i] = temp_holder[i];
    }
    for(int i = 0; i<dz; i++){
        nt_array[i] = temp_holder[dz+i];
    }
    for(int i = 0; i<dz; i++){
        np_array[i] = temp_holder[2*dz+i];
    }

    type_r[0] = RARE;
    type_r[dz-1] = UNSURR;
    for(int i = 1; i < dz-1; i++){
        type_r[i] = FIN;
    }

 
    








}