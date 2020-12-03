

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
#include "cmp.h"

#define PI 3.141592654

using namespace Lorene ;

/*
    This program generates the lapse function from the data produced by the python
    conversion code. The resulting LORENE scalar function is the same as would be
    saved in the GYOTO file. From about line 150 examples are given as to how to read
    the data from the scalar function given an r, phi and theta. 

*/



double * get_flatten_values(int size, char*);
double readdata(double *values, int l, int k, int j, int i, int nz, int nr, int nt, int np);
double readdata(double *values, int l, int k, int j, int i, int nz, int *nr, int *nt, int *np);
double *read_from_file(int size, char* name, int body, int it);

int main(int argc, char **argv) {

    /*
        Setup and reading of the files from the python conversion.
        Also shows a bit of how the LORENE objects are created
    */

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
    int nz = 7 ; 	// Number of domains
    /*
    int nr_array[]  = {55, 55, 55, 85, 17, 11};
    int nt_array[]  = {11, 11, 11, 11, 11, 11};
    int np_array[]  = {52, 72, 72, 82, 42, 42};
    */

    int nr_array[]  = {135, 135, 135, 135, 135, 67, 57};
    int nt_array[]  = {51, 51, 51, 51, 51, 51, 31};
    int np_array[]  = {142, 142, 142, 142, 122, 102, 62};

    // Gets the size of file written by the python program
    int size = 0;


    for(int j = 0; j < nz; j++){
      size += nr_array[j]*nt_array[j]*np_array[j];
    }



    int type_r[] = {RARE, FIN, FIN, FIN, FIN, FIN, UNSURR};
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = NONSYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified


    // Multi-domain grid construction:
    //Mg3d mgrid(nz, nr_array, type_r, nt_array, np_array, symmetry_theta, symmetry_phi, compact) ;
    Mg3d mgrid(nz, nr_array, type_r, nt_array, symmetry_theta, np_array, symmetry_phi, NULL) ;

    cout << mgrid << endl ;


    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    //double r_limits[] = {0., 0.5, 1.5, 4, 8, 20, __infinity} ;
    double r_limits[] = {0., 1., 1.5, 2.5 , 4, 8, 20, __infinity} ;

    Map_af map_1(mgrid, r_limits) ;   // Mapping construction
    Map_af map_2(mgrid, r_limits) ;   // Mapping construction

    // Sets the origin to one of the black holes
    
    map_1.set_ori(3, -0.0024, 0);
    map_2.set_ori(-3, 0.0024, 0);


    printf("%s\n", "[~] Reading Files");
    int it = 0;


    /* 
        Reads the laps functions from both of the black holes

    */
    // Lapse
    double *value_alp_1 = read_from_file(3 + size, "gxx", 1, it);
    double *value_alp_2 = read_from_file(3 + size, "gxx", 2, it);

    

    printf("%s\n", "[+] Files Successfully Read");

    //Allocation of quantities:

    // Lapse
    Scalar N_1(map_1) ;
    N_1.allocate_all() ;

    Scalar N_2(map_2) ;
    N_2.allocate_all() ;


    printf("%s\n", "[+] Quantities Successfully Allocated");
    
    
    // Assignes values from the read files
    for(int l=0;l<nz;l++){
        for(int k=0;k<np_array[l];k++){
            for(int j=0;j<nt_array[l];j++){
                for(int i=0;i<nr_array[l];i++){
                    N_1.set_grid_point(l, k, j, i) = readdata(value_alp_1, l, k, j, i,nz, nr_array, nt_array, np_array);
                    N_2.set_grid_point(l, k, j, i) = readdata(value_alp_2, l, k, j, i,nz, nr_array, nt_array, np_array);

                }
            }
        }
        }

    printf("%s\n", "[+] Quantities Successfully Filled");



    // Converts to spectral bases
    N_1.std_spectral_base();
    N_2.std_spectral_base();
    

    printf("%s\n", "[+] Quantities Successfully Made to Spectral Bases");


    /*
        #########################################
        The examples start here!
        #########################################
        
        First are some plots of the data. The lapse of the
        first black hole and the combined lapse are shown here.
    */


    double rmax=60;

    //des_coupe_z(N_1, 0., 4, "Lapse 1") ;
    //des_coupe_bin_z((Cmp)N_1, (Cmp)N_2, 0., -5, 5, -5, 5, "gxx",0x0,0x0,false) ;
    //des_coupe_bin_z((Cmp)N_1, (Cmp)N_2, 0., -20, 20, -20, 20, "gxx",0x0,0x0,false) ;
    /*
    des_coupe_z(N_2, 0., 4, "Lapse 2") ;
    des_coupe_bin_z((Cmp)N_1, (Cmp)N_2, 0., -10, 10, -10, 10, "Lapse",0x0,0x0,false) ;
    des_coupe_bin_z((Cmp)N_1, (Cmp)N_2, 0., -40, 40, -40, 40, "Lapse",0x0,0x0,false) ;
    arrete() ;
    */




    /*
        This gets the position of the first black hole in cartesian coordinates from the map object
    */

    double ori_x = map_1.get_ori_x();
    double ori_y = map_1.get_ori_y();
    double ori_z = map_1.get_ori_z();
    
    cout << "Position of black hole 1 is (" << ori_x << ", " << ori_y << ", " << ori_z << ")" << endl;

    /*
        The value of the lapse function can be retrived for any point r, theta, phi.
        Sadly the orgin of this coordinates system is at the black holes (the origin of the map object).
    */


    // Lapse function at the center of the BHs
    cout << N_1.val_point(0,0,0) << endl;
    cout << N_2.val_point(0,0,0) << endl;
    

    /*
        A small loop to print out the value of the lapse of the firts BH.
        This crosses the second BH and the origin, and if plotted shows
        the profile and the effect of the smooting function.
    */
    double N1, N2, r, phi;
    
    for(double x = -10; x <= 10; x += 0.005){
        r = sqrt(x*x );
        phi = atan2(x,0) + PI;
        cout  << x << "   " << N_2.val_point(r,0,phi) << endl;
    }


    
    
    return EXIT_SUCCESS ;


}



/* 
    ############################################################
    Functions for reading the files from the python conversion
    Can be ignored!
    ############################################################
*/



double readdata(double *values, int l, int k, int j, int i, int nz, int nr, int nt, int np){
    int index = 3 + l*nr*nt*np + k*nt*nr + j*nr + i;
    return *(values + index);
}

double readdata(double *values, int l, int k, int j, int i, int nz, int *nr, int *nt, int *np){
    int outer_size = 0;

    for(int m = 0; m < l; m++){
      outer_size += nr[m]*nt[m]*np[m];
    }

    int index = 3 + outer_size + k*nt[l]*nr[l] + j*nr[l] + i;

    //cout << index << endl;
    return *(values + index);
}

double *get_flatten_values(int size, char* filename){

    const int max_size = 2048*128*64;

    FILE *myFile;
    //myFile = fopen("flatten.txt", "r");
    myFile = fopen(filename, "r");

    if(size > max_size){
        cout << "Too large a size for the flatten data size!" << endl;
        exit(1);
    }

    /*
    //read file into array
    static double values[max_size];
    int i;
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (1);
    }

    for (i = 0; i < size; i++){
        fscanf(myFile, "%lf,", &values[i] );
    }

    */


    double* values = static_cast<double *>(malloc(sizeof(double) * size));
    int i;
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (1);
    }

    for (i = 0; i < size; i++){
        fscanf(myFile, "%lf,", &values[i] );
    }
    
    fclose(myFile);
    /*
    for (i = 0; i < size; i++){
      printf("Number is: %f\n\n", values[i]);
    }
    */

    return values;

}

double *read_from_file(int size, char* name, int body, int it){
    char file[20];

    sprintf(file, "%s_%d_body%d.txt", name, it, body);
    double *value = get_flatten_values(size, file);

    return value;
}
