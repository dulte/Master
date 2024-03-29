

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


using namespace Lorene ;
double * get_flatten_values(int size, char*);
double readdata(double *values, int l, int k, int j, int i, int nz, int nr, int nt, int np);
double readdata(double *values, int l, int k, int j, int i, int nz, int *nr, int *nt, int *np);
double *read_from_file(int size, char* name, int body, int it);

int main(int argc, char **argv) {

    if(argc < 5){
        cout << "Usage: program x_origin y_origin z_origin 0 (write to screen)/1 (read values from file) body iteration" << endl;
        exit(1);
    }


    /*
    #######################################################
        User defined variables  
        Change to change conversion!
    #######################################################
    */
    int nz = 6 ; 	// Number of domains

    // Domain Resolutions
    int nr_array[]  = {25, 25, 25, 25, 25, 25};
    int nt_array[]  = {7,7,7, 7,7,7};
    int np_array[]  = {4,4,4, 4,4,4};

    //Type of Domain
    int type_r[] = {RARE, FIN, FIN,FIN,FIN, UNSURR};

    // Domain Limits
    double r_limits[] = {0.,0.51, 1, 2, 4, 8, __infinity} ;

    /*
    #######################################################
        Rest of code...
    #######################################################
    */




    int size = 0;
    for(int j = 0; j < nz; j++){
      size += nr_array[j]*nt_array[j]*np_array[j];
    }



    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = NONSYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified


    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr_array, type_r, nt_array, symmetry_theta, np_array, symmetry_phi, NULL) ;

    cout << mgrid << endl ;









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

    /*
        This part is for getting the collocation points
    */

    if(read_write == 0){
        cout << "+" << endl;
        cout << x << endl;

        cout << "+" << endl;
        cout << y << endl;

        cout << "+" << endl;
        cout << z << endl;









    /*
        This part is for doing the spectral transformation
    */


    }else if(read_write == 1){
        int body;
        
        sscanf(argv[5], "%d", &body);

        int it;
        sscanf(argv[6], "%d", &it);

        /* This is horrble practice, and I should have made just one file.
        To be fixed...
        */

        printf("%s\n", "[~] Reading Files");

        // Lapse
        double *value_alp = read_from_file(3 + size, "alp", body, it);

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


        printf("%s\n", "[+] Files Successfully Read");

        //Allocation of quantities:
        
        // Lapse
        Scalar N(map) ;
        N.allocate_all() ;

        // Shift
        Vector beta(map,CON,map.get_bvect_cart());
        beta.allocate_all();

        // Metric
        Sym_tensor gamma(map,COV,map.get_bvect_cart());
        gamma.allocate_all();

        // Curavture
        Sym_tensor curvature(map,COV,map.get_bvect_cart());
        curvature.allocate_all();

        printf("%s\n", "[+] Quantities Successfully Allocated");



        // Assignes values from the read files

        for(int l=0;l<nz;l++){
            for(int k=0;k<np_array[l];k++){
                for(int j=0;j<nt_array[l];j++){
                    for(int i=0;i<nr_array[l];i++){
                        //cout << l << " " << k << " " << j << " " << i << " " << l*k*j*i << " " << size << endl;
                        double Ndata = readdata(value_alp, l, k, j, i,nz, nr_array, nt_array, np_array);
                        N.set_grid_point(l, k, j, i) = Ndata;
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
                    }
                }
            }
    	  }

        printf("%s\n", "[+] Quantities Successfully Filled");



        
        



        // Converts to spectral bases
        N.std_spectral_base();
        beta.std_spectral_base();
        gamma.std_spectral_base();
        curvature.std_spectral_base();



        printf("%s\n", "[+] Quantities Successfully Made to Spectral Bases");

        Sym_tensor cart_gamma(map,COV,map.get_bvect_cart());
        cart_gamma = gamma;

        

        // Converts from cartesian to spherical
        
        beta.change_triad(map.get_bvect_spher());
        gamma.change_triad(map.get_bvect_spher());
        curvature.change_triad(map.get_bvect_spher());
        

        printf("%s\n", "[+] Quantities Successfully Converted");
        

        
        Metric metric(map.flat_met_spher());
        metric = gamma;


        // Plotting for testing
        /*
        double rmax=8;
        des_meridian(gamma(1,1), 0.51, rmax, "g_{xx}", 3) ;
        des_meridian(N, 0.5, rmax, "\\alpha", 4) ;
        des_meridian(gamma(1,1) - pow((1+1/(2*map.r)),4), 0.52, rmax, "Simulated - Analytical g_{xx}", 5) ;
        des_meridian(N - (1-1/(2*map.r))/(1+1/(2*map.r)), 0.52, rmax, "Simulated - Analytical \\alpha", 6) ;
        
        des_coupe_z(gamma(1,1), 1, nz-1, "g_{xx}") ;
        des_coupe_z(N, 1, nz-1, "\\alpha") ;

        arrete() ;
        */

        // Saves to Gyoto Readable File
        
        char out_name[20];
        sprintf(out_name, "bbh_%d_body%d.d",it, body);


        FILE* file_out = fopen(out_name, "w") ;
        double total_time = 0. ; // for compatibility

        fwrite_be(&total_time, sizeof(double), 1, file_out) ;
        mgrid.sauve(file_out) ;
        map.sauve(file_out) ;
        N.sauve(file_out) ;
        beta.sauve(file_out) ;
        Metric(gamma).cov().sauve(file_out) ;
        Metric(gamma).con().sauve(file_out);
        
        curvature.sauve(file_out) ;

        fclose(file_out) ;

        printf("%s\n", "[+] Quantities Successfully Written to File");

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

double readdata(double *values, int l, int k, int j, int i, int nz, int *nr, int *nt, int *np){
    int outer_size = 0;

    for(int m = 0; m < l; m++){
      outer_size += nr[m]*nt[m]*np[m];
    }

    int index = 3 + outer_size + k*nt[l]*nr[l] + j*nr[l] + i;

    return *(values + index);
}

double *get_flatten_values(int size, char* filename){

    const int max_size = 2048*128*64;

    FILE *myFile;
    myFile = fopen(filename, "r");

    if(size > max_size){
        cout << "Too large a size for the flatten data size!" << endl;
        exit(1);
    }



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

    return values;

}

double *read_from_file(int size, char* name, int body, int it){
    char file[20];

    sprintf(file, "%s_%d_body%d.txt", name, it, body);
    double *value = get_flatten_values(size, file);

    return value;
}
