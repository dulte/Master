

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
double *read_from_file(int size, char* name, int body, int it);

int main(int argc, char **argv) {

    if(argc < 6){
        cout << "Usage: program x_origin y_origin z_origin 0 (write to screen)/1 (read values from file) body iteration" << endl;
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


    //sscanf(argv[5], "%s", filename);

    if(read_write == 0){
        cout << "+" << endl;
        cout << x << endl;

        cout << "+" << endl;
        cout << y << endl;

        cout << "+" << endl;
        cout << z << endl;
    }else if(read_write == 1){
        int body;
        sscanf(argv[5], "%d", &body);

        int it;
        sscanf(argv[6], "%d", &it);

        /* This is horrble practis, and I should have made just one file.
        To be fixed...
        */

        // Lapse
        double *value_alp = read_from_file(3 + nz*nr*nt*np, "alp", body, it);

        // Shift
        double *value_betax = read_from_file(3 + nz*nr*nt*np, "betax", body, it);
        double *value_betay = read_from_file(3 + nz*nr*nt*np, "betay", body, it);
        double *value_betaz = read_from_file(3 + nz*nr*nt*np, "betaz", body, it);

        // Metric
        double *value_gxx = read_from_file(3 + nz*nr*nt*np, "gxx", body, it);
        double *value_gxy = read_from_file(3 + nz*nr*nt*np, "gxy", body, it);
        double *value_gxz = read_from_file(3 + nz*nr*nt*np, "gxz", body, it);

        double *value_gyy = read_from_file(3 + nz*nr*nt*np, "gyy", body, it);
        double *value_gyz = read_from_file(3 + nz*nr*nt*np, "gyz", body, it);
        double *value_gzz = read_from_file(3 + nz*nr*nt*np, "gzz", body, it);

        // Curvature
        double *value_kxx = read_from_file(3 + nz*nr*nt*np, "kxx", body, it);
        double *value_kxy = read_from_file(3 + nz*nr*nt*np, "kxy", body, it);
        double *value_kxz = read_from_file(3 + nz*nr*nt*np, "kxz", body, it);

        double *value_kyy = read_from_file(3 + nz*nr*nt*np, "kyy", body, it);
        double *value_kyz = read_from_file(3 + nz*nr*nt*np, "kyz", body, it);
        double *value_kzz = read_from_file(3 + nz*nr*nt*np, "kzz", body, it);





        // Lapse
        Scalar N(map) ;
        N.allocate_all() ;

        // Shift

        // Metric

        // Curavture



        for(int l=0;l<nz;l++){
            for(int k=0;k<np;k++){
                for(int j=0;j<nt;j++){
                    for(int i=0;i<nr;i++){
                        double Ndata = readdata(value_alp, l, k, j, i, nz, nr, nt, np);
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

double *read_from_file(int size, char* name, int body, int it){
    char file[20];
    sprintf(file, "%s_%d_body%d.txt", name, body,it);
    double *value = get_flatten_values(size, file);

    return value;
}
