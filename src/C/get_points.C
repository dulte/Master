

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
double * get_flatten_values(int size, char*);
double readdata(double *values, int l, int k, int j, int i, int nz, int nr, int nt, int np);

int main(int argc, char **argv) {

    if(argc < 6){
        cout << "Usage: program x_origin y_origin z_origin 0 (write to screen)/1 (read values from file) filename" << endl;
        exit(1);
    }

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------

    int nz = 3 ; 	// Number of domains
    int nr = 7 ; 	// Number of collocation points in r in each domain
    int nt = 5 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified

    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

    cout << mgrid << endl ;


    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 1., 2., __infinity} ;

    Map_af map(mgrid, r_limits) ;   // Mapping construction

    // Sets the origin to one of the black holes
    double xo, yo, zo;
    sscanf(argv[1], "%lf", &xo);
    sscanf(argv[2], "%lf", &yo);
    sscanf(argv[3], "%lf", &zo);
    map.set_ori(xo, yo, zo);

    cout << map << endl ;

    // Denomination of various coordinates associated with the mapping
    // ---------------------------------------------------------------

    const Coord& x = map.xa ;        // x field
    const Coord& y = map.ya ;        // y field
    const Coord& z = map.za ;        // z field

    int read_write;
    sscanf(argv[4], "%d", &read_write);

    char* filename = argv[5];
    //sscanf(argv[5], "%s", filename);

    if(read_write == 0){
        cout << "+" << endl;
        cout << x << endl;

        cout << "+" << endl;
        cout << y << endl;

        cout << "+" << endl;
        cout << z << endl;
    }else if(read_write == 1){
        double *value_array = get_flatten_values(3 + nz*nr*nt*np, filename);
        cout << readdata(value_array, 0, 1, 3, 4, nz, nr, nt, np) << endl;

        //Lapse
        Scalar N(map) ;
        N.allocate_all() ;

        for(int l=0;l<nz;l++){
            for(int k=0;k<np;k++){
                for(int j=0;j<nt;j++){
                    for(int i=0;i<nr;i++){
                        double Ndata = readdata(value_array, l, k, j, i, nz, nr, nt, np);
                        N.set_grid_point(l, k, j, i) = Ndata;
                    }
                }
            }
    	}
        N.std_spectral_base();
        cout << "hello" << endl;
        double rmax=4;
        //des_meridian(N, 0, rmax, "N", 1) ;



        //des_coupe_z(N, 0., 2, "Lapse") ;


        //arrete() ;

        cout << "[+] Successfully Read From File and Saved Quantity" << endl;
    }

    else{
        cout << "Fourth argument needs to be 1 or 0!" << endl;
        exit(1);
    }




    return EXIT_SUCCESS ;
}

double readdata(double *values, int l, int k, int j, int i, int nz, int nr, int nt, int np){
    int index = 3 + l*nr*nt*np + k*nt*nr + j*nr + i;
    return *(values + index);
}


double *get_flatten_values(int size, char* filename){


    FILE *myFile;
    //myFile = fopen("flatten.txt", "r");
    myFile = fopen(filename, "r");

    if(size > 1024){
        cout << "Too large a size for the flatten data size!" << endl;
        exit(1);
    }

    //read file into array
    static double values[1024];
    int i;
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (1);
    }

    for (i = 0; i < size; i++){
        fscanf(myFile, "%lf,", &values[i] );
    }

    /*
    for (i = 0; i < size; i++){
        printf("Number is: %f\n\n", values[i]);
    }
    */


    fclose(myFile);

    return values;


}
