

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


/*
TODO:
    - Use get_bvect_cart() since this is read as cartesian!
      Then use change_triad(map.get_bvect_spher())

*/



using namespace Lorene ;
double * get_flatten_values(int size, char*);
double readdata(double *values, int l, int k, int j, int i, int nz, int nr, int nt, int np);
double readdata(double *values, int l, int k, int j, int i, int nz, int *nr, int *nt, int *np);
double *read_from_file(int size, char* name, int body, int it);

int main(int argc, char **argv) {

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
    int nz = 6 ; 	// Number of domains
    /*
    int nr = 25; 	// Number of collocation points in r in each domain
    int nt = 11 ; 	// Number of collocation points in theta in each domain
    int np = 42 ; 	// Number of collocation points in phi in each domain
    */

    int nr_array[]  = {135, 135, 135, 135, 67, 57};
    int nt_array[]  = {51, 51, 51, 51, 51, 31};
    int np_array[]  = {142, 142, 142, 122, 102, 62};



    // int size = nz*nr*np*nt
    int size = 0;


    for(int j = 0; j < nz; j++){
      size += nr_array[j]*nt_array[j]*np_array[j];
    }



    int type_r[] = {RARE, FIN, FIN, FIN, FIN, UNSURR};
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
    double r_limits[] = {0., 0.5, 1.5, 4, 8, 20, __infinity} ;

    Map_af map_1(mgrid, r_limits) ;   // Mapping construction
    Map_af map_2(mgrid, r_limits) ;   // Mapping construction

    // Sets the origin to one of the black holes
    
    map_1.set_ori(3, -0.0024, 0);
    map_2.set_ori(-3, 0.0024, 0);


    /*
    cout << map << endl ;

    // Denomination of various coordinates associated with the mapping
    // ---------------------------------------------------------------

    const Coord& x = map.xa ;        // x field
    const Coord& y = map.ya ;        // y field
    const Coord& z = map.za ;        // z field


    */


    //sscanf(argv[5], "%s", filename);

    
    
    

    int it = 0;
    

    /* This is horrble practice, and I should have made just one file.
    To be fixed...
    */

    printf("%s\n", "[~] Reading Files");

    // Lapse
    double *value_alp_1 = read_from_file(3 + size, "alp", 1, it);
    double *value_alp_2 = read_from_file(3 + size, "alp", 2, it);

    /*
    // Shift
    double *value_betax = read_from_file(3 + size, "betax", body, it);
    double *value_betaz = read_from_file(3 + size, "betaz", body, it);
    double *value_betay = read_from_file(3 + size, "betay", body, it);

    // Metric
    double *value_gxx = read_from_file(3 + size, "gxx", body, it);
    double *value_gxy = read_from_file(3 + size, "gxy", body, it);
    double *value_gxz = read_from_file(3 + size, "gxz", body, it);

    double *value_gyy = read_from_file(3 + size, "gyy", body, it);
    double *value_gyz = read_from_file(3 + size, "gyz", body, it);
    double *value_gzz = read_from_file(3 + size, "gzz", body, it);

    // Curvature
    double *value_kxx = read_from_file(3 + size, "kxx", body, it);
    double *value_kxy = read_from_file(3 + size, "kxy", body, it);
    double *value_kxz = read_from_file(3 + size, "kxz", body, it);

    double *value_kyy = read_from_file(3 + size, "kyy", body, it);
    double *value_kyz = read_from_file(3 + size, "kyz", body, it);
    double *value_kzz = read_from_file(3 + size, "kzz", body, it);
    */

    printf("%s\n", "[+] Files Successfully Read");

    //Allocation of quantities:

    // Lapse
    Scalar N_1(map_1) ;
    N_1.allocate_all() ;

    Scalar N_2(map_2) ;
    N_2.allocate_all() ;

    /*

    // Shift
    Vector beta(map,CON,map.get_bvect_cart());
    beta.allocate_all();

    // Metric
    Sym_tensor gamma(map,COV,map.get_bvect_cart());
    gamma.allocate_all();

    // Curavture
    Sym_tensor curvature(map,COV,map.get_bvect_cart());
    curvature.allocate_all();
    */

    printf("%s\n", "[+] Quantities Successfully Allocated");
    // Assignes values from the read files


    for(int l=0;l<nz;l++){
        for(int k=0;k<np_array[l];k++){
            for(int j=0;j<nt_array[l];j++){
                for(int i=0;i<nr_array[l];i++){
                    //cout << l << " " << k << " " << j << " " << i << " " << l*k*j*i << " " << size << endl;
                    //double Ndata = readdata(value_alp, l, k, j, i,nz, nr_array, nt_array, np_array);
                    N_1.set_grid_point(l, k, j, i) = readdata(value_alp_1, l, k, j, i,nz, nr_array, nt_array, np_array);
                    N_2.set_grid_point(l, k, j, i) = readdata(value_alp_2, l, k, j, i,nz, nr_array, nt_array, np_array);

                    /*
                    beta.set(1).set_grid_point(l, k, j, i) = readdata(value_betax, l, k, j, i, nz, nr_array, nt_array, np_array);
                    beta.set(2).set_grid_point(l, k, j, i) = readdata(value_betay, l, k, j, i, nz, nr_array, nt_array, np_array);
                    beta.set(3).set_grid_point(l, k, j, i) = readdata(value_betaz, l, k, j, i, nz, nr_array, nt_array, np_array);

                    gamma.set(1,1).set_grid_point(l, k, j, i) = readdata(value_gxx, l, k, j, i, nz, nr_array, nt_array, np_array);
                    gamma.set(1,2).set_grid_point(l, k, j, i) = readdata(value_gxy, l, k, j, i, nz, nr_array, nt_array, np_array);
                    gamma.set(1,3).set_grid_point(l, k, j, i) = readdata(value_gxz, l, k, j, i, nz, nr_array, nt_array, np_array);

                    gamma.set(2,2).set_grid_point(l, k, j, i) = readdata(value_gyy, l, k, j, i, nz, nr_array, nt_array, np_array);
                    gamma.set(2,3).set_grid_point(l, k, j, i) = readdata(value_gyz, l, k, j, i, nz, nr_array, nt_array, np_array);
                    gamma.set(3,3).set_grid_point(l, k, j, i) = readdata(value_gzz, l, k, j, i, nz, nr_array, nt_array, np_array);

                    gamma.set(2,1).set_grid_point(l, k, j, i) = readdata(value_gxy, l, k, j, i, nz, nr_array, nt_array, np_array);
                    gamma.set(3,1).set_grid_point(l, k, j, i) = readdata(value_gxz, l, k, j, i, nz, nr_array, nt_array, np_array);
                    gamma.set(3,2).set_grid_point(l, k, j, i) = readdata(value_gyz, l, k, j, i, nz, nr_array, nt_array, np_array);


                    curvature.set(1,1).set_grid_point(l, k, j, i) = readdata(value_kxx, l, k, j, i, nz, nr_array, nt_array, np_array);
                    curvature.set(1,2).set_grid_point(l, k, j, i) = readdata(value_kxy, l, k, j, i, nz, nr_array, nt_array, np_array);
                    curvature.set(1,3).set_grid_point(l, k, j, i) = readdata(value_kxz, l, k, j, i, nz, nr_array, nt_array, np_array);

                    curvature.set(2,2).set_grid_point(l, k, j, i) = readdata(value_kyy, l, k, j, i, nz, nr_array, nt_array, np_array);
                    curvature.set(2,3).set_grid_point(l, k, j, i) = readdata(value_kyz, l, k, j, i, nz, nr_array, nt_array, np_array);
                    curvature.set(3,3).set_grid_point(l, k, j, i) = readdata(value_kzz, l, k, j, i, nz, nr_array, nt_array, np_array);

                    curvature.set(2,1).set_grid_point(l, k, j, i) = readdata(value_kxy, l, k, j, i, nz, nr_array, nt_array, np_array);
                    curvature.set(3,1).set_grid_point(l, k, j, i) = readdata(value_kxz, l, k, j, i, nz, nr_array, nt_array, np_array);
                    curvature.set(3,2).set_grid_point(l, k, j, i) = readdata(value_kyz, l, k, j, i, nz, nr_array, nt_array, np_array);
                    */
                }
            }
        }
        }

    printf("%s\n", "[+] Quantities Successfully Filled");








    // Converts to spectral bases
    N_1.std_spectral_base();
    N_2.std_spectral_base();
    //beta.std_spectral_base();
    //gamma.std_spectral_base();
    //curvature.std_spectral_base();
    


    printf("%s\n", "[+] Quantities Successfully Made to Spectral Bases");
    double rmax=60;

    //des_meridian(N_1, 0, rmax, "N", 1) ;
    //des_meridian(N_1, 0, 20, "N", 1) ;
    des_coupe_z(N_1, 0., 4, "Lapse 1") ;
    des_coupe_z(N_2, 0., 4, "Lapse 2") ;
    des_coupe_bin_z((Cmp)N_1, (Cmp)N_2, 0., -20, 20, -20, 20, "Lapse",0x0,0x0,false) ;
    des_coupe_bin_z((Cmp)N_1, (Cmp)N_2, 0., -10, 10, -10, 10, "Lapse",0x0,0x0,false) ;
    des_coupe_bin_z((Cmp)N_1, (Cmp)N_2, 0., -40, 40, -40, 40, "Lapse",0x0,0x0,false) ;
    arrete() ;






    return EXIT_SUCCESS ;
    // Converts from cartesian to spherical
    /*
    cout << map.get_bvect_spher() << endl;
    beta.change_triad(map.get_bvect_spher());
    gamma.change_triad(map.get_bvect_spher());
    curvature.change_triad(map.get_bvect_spher());
    

    printf("%s\n", "[+] Quantities Successfully Converted");
    */


    // Plotting for testing
    cout << N_1.val_point(3,0,1) << endl;

    FILE *fptr;
    double N1, N2, r, phi;

    if ((fptr = fopen("N_spherical.txt","w")) == NULL){
        printf("Error! opening file");

        // Program exits if the file pointer returns NULL.
        exit(1);
    }
    
    for(double x = -7; x<=7; x += 0.02)
    {
        cout << x << endl;
        for(double y = -7; y<=7; y += 0.02){   
            r = sqrt((x+3)*(x+3) + y*y);// + (y-0.0024)*(y-0.0024));
            //phi = atan2(y-0.0024,x+3);
            phi = atan2(x+3,y) + PI;
            N1 = N_1.val_point(r,phi,0); //+ N_2.val_point(r,phi,0);
            r = sqrt((x-3)*(x-3) + y*y);// + (y+0.0024)*(y+0.0024));
            //phi = atan2(y+0.0024,x-3);
            phi = atan2(x-3,y) + PI;

            N2 = N_2.val_point(r,phi,0);

            //fwrite(&N, sizeof(double), 1, fptr); 
            fprintf(fptr, "%f %f %f %f\n", N1, N2, x, y);
        }
    }
    fclose(fptr); 

    
    //des_meridian(N_2, 0, rmax, "N", 1) ;
    des_coupe_z(N_1, 0., 4, "Lapse") ;
    des_coupe_z(N_2, 0., 4, "Lapse") ;

    des_coupe_z(N_1+N_2, 0., 4, "Lapse") ;
    //des_coupe_z(gamma(1,1), 0., 4, "g_{xx}") ;
    //des_coupe_z(curvature(1,1), 0., 5, "K_{xx}") ;
    arrete() ;


    // Saves to Gyoto Readable File
    /*Make correct file name... */



    printf("%s\n", "[+] Quantities Successfully Written to File");

    cout << "[+] Successfully Read From File and Saved Quantity" << endl;
    




    return EXIT_SUCCESS ;
}

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
