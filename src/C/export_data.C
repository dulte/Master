/*
 *  Export data to Lorene in a format readable by Gyoto
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include "stdlib.h"
#include "assert.h"
#include "math.h"

// Lorene headers
#include "nbr_spx.h"
#include "tensor.h"
#include "metric.h"
#include "graphique.h"
#include "utilitaires.h"

#define Rmax 300
#define rh 15.72646372
#define m 0.06

using namespace Lorene ;

void init(double f[][3], double f2[][2], const int n){
    double p,sig;
    double u[n][2];

    for(int i=0;i<2;i++){
        f2[0][i]=0.0;
        u[0][i]=0.0;
    }

    for(int j=0;j<2;j++){
        for(int i=1;i<n-1;i++){
            sig=(f[i][0]-f[i-1][0])/(f[i+1][0]-f[i-1][0]);
            p=sig*f2[i-1][j]+2.0;
            f2[i][j]=(sig-1.0)/p;
            u[i][j]=(f[i+1][j+1]-f[i][j+1])/(f[i+1][0]-f[i][0])-(f[i][j+1]-f[i-1][j+1])/(f[i][0]-f[i-1][0]);
            u[i][j]=(6.0*u[i][j]/(f[i+1][0]-f[i-1][0])-sig*u[i-1][j])/p;
        }
        f2[n-1][j]=0.0;
        for(int k=n-2;k>=0;k--)   f2[k][j]=f2[k][j]*f2[k+1][j]+u[k][j];
    }
}

double inter(double f[][3], double f2[][2], const int n, double r, int c){
    int k;
    double h,b,a;

    int klo=0;
    int khi=n-1;
    while(khi-klo > 1){
        k=(khi+klo)>>1;
        if(f[k][0]>r)   khi=k;
        else klo=k;
    }
    h=f[khi][0]-f[klo][0];
    a=(f[khi][0]-r)/h;
    b=(r-f[klo][0])/h;
    return a*f[klo][c]+b*f[khi][c]+((a*a*a-a)*f2[klo][c-1]+(b*b*b-b)*f2[khi][c-1])*(h*h)/6.0;
}

int main() {

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------

    int nz = 5 ; 	// Number of domains
    int nr = 33 ; 	// Number of collocation points in r in each domain
    int nt = 1 ; 	// Number of collocation points in theta in each domain
    int np = 1 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified

    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

    cout << mgrid << endl ;


    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., /*rh+0.1*/6.1, /*rh+1.0*/7., /*32.*/12., /*64.*/24., __infinity} ;

    Map_af map(mgrid, r_limits) ;   // Mapping construction

    cout << map << endl ;

    // Denomination of various coordinates associated with the mapping
    // ---------------------------------------------------------------

    const Coord& xx = map.x ;        // r field
    const Coord& yy = map.y ;        // r field
    const Coord& zz = map.z ;        // r field

    // to be called as xx(l, k, j, i)
    cout << rr << endl;

    //Lecture du fichier de données
    int nblinestemp=0;
    string line;
    double temp;
    ifstream bh("bh.out");
    ifstream par("par.out");

    while(getline(bh,line)) ++nblinestemp;
    const int n=nblinestemp;

    bh.clear();
    bh.seekg(0,ifstream::beg);
    double f[n][3];
    double f2[n][2];

    for(int i=0;i<n;i++){
        for(int j=0;j<3;j++){
            bh >> f[i][j];
            if(j<2) f2[i][j]=0.0;
        }
    }

    init(f,f2,n);
    double r=6.00001;
    double A,B;
    par >> A;
    par >> B;

    //Affichage interpolation
    /*ofstream out("inter.out",ios::out | ios::trunc);
    while(r<Rmax){
        cout << r << endl;
        out << r << "\t" << 1/pow(inter(f,f2,n,r,1),2) << "\t" << inter(f,f2,n,r,2) << endl;
        r+=0.01;
    }*/

    //Lapse
    Scalar N(map) ;
    N.allocate_all() ;
    N.set_domain(0)=1 ; // set to zero in first domain
    for(int l=1;l<nz;l++){
    	for(int k=0;k<np;k++){
    	    for(int j=0;j<nt;j++){
    	        for(int i=0;i<nr;i++){
            xg = (*(xx.c))(l, k, j, i)
            yg = yy(l, k, j, i)
            zg = zz(l, k, j, i)

           Ndata = readdata(xg, yg, zg, ...)
    		N.set_grid_point(l, k, j, i) = Ndata
    		}
    	}
    }
    N.std_spectral_base();

    //Shift vector
    //Use get_bvect_cart() since this is read as cartesian!
    // Then use change_triad(map.get_bvect_spher())
    Vector beta(map,CON,map.get_bvect_spher());
    beta.set_etat_zero();

    //3-metric
    Sym_tensor gamma(map,COV,map.get_bvect_spher());
    gamma.allocate_all();
    gamma.set(1,1).set_domain(0)=1; // set to zero in first domain
    for(int i=1;i<nz;i++){
    	for(int j=0;j<nr;j++){
    		r=(*(rr.c))(i,0,0,j);
    		if(r < Rmax) gamma.set(1,1).set_grid_point(i,0,0,j)=1/pow(inter(f,f2,n,r,1)/*sqrt(1.0-rh/r)*/,2);
    		else{
    			if(i==nz-1 && j==nr-1) gamma.set(1,1).set_grid_point(i,0,0,j)=1.0;
    			else gamma.set(1,1).set_grid_point(i,0,0,j)=1/pow(1.0-(A/(2.0*r))+B*(m*r+1.0)*exp(-m*r)/(2.0*r),2);
    		}
    	}
    }
    gamma.set(1,2)=0;
    gamma.set(1,3)=0;
    gamma.set(2,2)=1;
    gamma.set(2,3)=0;
    gamma.set(3,3)=1;
    gamma.std_spectral_base();

    gam =Metric(gamma)
    gam.con()

    Sym_tensor inv_gamma(map,CON,map.get_bvect_spher());
    inv_gamma.allocate_all();
    inv_gamma.set(1,1).set_domain(0)=1; // set to zero in first domain
    for(int i=1;i<nz;i++){
    	for(int j=0;j<nr;j++){
    		r=(*(rr.c))(i,0,0,j);
    		if(r < Rmax) inv_gamma.set(1,1).set_grid_point(i,0,0,j)=pow(inter(f,f2,n,r,1)/*sqrt(1.0-rh/r)*/,2);
    		else{
    			if(i==nz-1 && j==nr-1) inv_gamma.set(1,1).set_grid_point(i,0,0,j)=1.0;
    			else inv_gamma.set(1,1).set_grid_point(i,0,0,j)=pow(1.0-(A/(2.0*r))+B*(m*r+1.0)*exp(-m*r)/(2.0*r),2);
    		}
    	}
    }
    inv_gamma.set(1,2)=0;
    inv_gamma.set(1,3)=0;
    inv_gamma.set(2,2)=1;
    inv_gamma.set(2,3)=0;
    inv_gamma.set(3,3)=1;
    inv_gamma.std_spectral_base();

    //Extrinsic curvature
    Sym_tensor kk(map,COV,map.get_bvect_spher());
    kk.set_etat_zero();

    //Angular momentum
    double aa=0;

    //Drawings
    double rmax=1000;
    des_meridian(N, 0, rmax, "N", 1) ;
    des_meridian(gamma(1,1), 0, rmax, "gamma_11", 4) ;

	des_coef_xi(N.get_spectral_va(),1,0,0);
	des_coef_xi(N.get_spectral_va(),2,0,0);
	des_coef_xi(gamma(1,1).get_spectral_va(),1,0,0);
	des_coef_xi(gamma(1,1).get_spectral_va(),2,0,0);

    arrete() ;

    //Output file
    FILE* file_out = fopen("gyoto_hairy_bh.d", "w") ;
    double total_time = 0. ; // for compatibility

    fwrite_be(&total_time, sizeof(double), 1, file_out) ;
    mgrid.sauve(file_out) ;
    map.sauve(file_out) ;
    N.sauve(file_out) ;
    beta.sauve(file_out) ;
    gamma.sauve(file_out) ;
    inv_gamma.sauve(file_out) ;
    gam.con().sauve
    kk.sauve(file_out) ;
    fwrite_be(&aa, sizeof(double), 1, file_out) ;

    fclose(file_out) ;

    return EXIT_SUCCESS ;
}
