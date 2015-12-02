/********************************************************************/
/*	Program to calculate constant-temperature level densities of    */
/*	the format used in the TALYS tabulated level densities,         */
/*	ldmodel 4, for 89Y                                              */
/*																	*/
/*	Version Sunday 8 February 2015, Cecilie                         */
/********************************************************************/

/* 89Y */
/* LOWER LIMIT, D0_LOWER = 142.65 eV, FG09 SPIN CUTOFF E&B 2009 */

using namespace std;

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

/* Declare parameters */
const double T_ct  = 1.050000;      // temperature parameter
const double E0    = 0.097726;  // shift for the CT level density
 
//const double a_par = 8.856;     // level density parameter
const double A     = 89;      // mass A
//const double E1    = 0.277;     // Fermi-gas shift for the spin cutoff parameter, E&B 2005/2006
const double Pa = 1.316; //  deuteron pairing shift Pa' for the spin cutoff parameter, E&B 2009


/* Declare functions for constant-temperature level density and spin cutoff parameter */
double ct_model(double );
double spin_distribution(double , double );

int main(){
    
    /* Open  input and output files */
	ifstream inputEx89Yfile("ex_and_T_89Y.txt");
	FILE *output89Yfile;
    output89Yfile = fopen("ct_89Y_ldmodel4_format.txt","w");
    
	/* Vectors for calculated level densities. */
    double e_calc[56], temp_calc[56], nld_calc[56];
    double start_spin = 0.5;    // the first spin, 1/2 (odd) or 0 (even). Loop up to 59/2 = 29.5 if odd spin, up to 29 (starting from 0) if even spin
    int stop_spin = 30; // spin loop stops after 30 iterations
    double n_cum = 0.;  // cumulative number of levels. Remember different bin size on Ex!
    double ex_bin1 = 0.25; // 0.25 MeV from Ex=0.25 - 5.00 MeV, i=0-19
    double ex_bin2 = 0.50; // 0.50 MeV from Ex=5.50 - 10.0 MeV, i=20-29
    double ex_bin3 = 1.00; // 1.00 MeV from Ex=11.0 - 20.0 MeV, i=30-39
    double ex_bin4 = 2.50; // 2.50 MeV from Ex=22.5 - 25.0 MeV, i=40-41
    double ex_bin5 = 5.00; // 5.00 MeV from Ex=25.0 - 30.0 MeV, i=41-42
    double ex_bin6 = 10.0; // 10.0 MeV from Ex=30.0 - 150. MeV, i=43-54
    // Dummy variables
    double Ex;
    double I;
    double levdens;
    
    string line;
    line.resize(300);
    
	/* Read excitation energy and temperature from the ldmodel 4 file */
    int i=0;
    double x, y, z;
    
    
    getline(inputEx89Yfile,line);
    getline(inputEx89Yfile,line);
    i=0;
    while(inputEx89Yfile){
        inputEx89Yfile >> x >> y;
        e_calc[i]    = x;
        temp_calc[i] = y;
        //cout << e_calc[i] << " " << temp_calc[i] << endl;
        i++;
    }
     
    /* Calculate ct level density and print to file */
    /* For each Ex, calculate total level density, cumulative level density, and spin-dependent level density */
	for(i=0; i<55; i++){
        nld_calc[i] = ct_model(e_calc[i]);

        if(i<20)
            n_cum += nld_calc[i]*ex_bin1;
        if(i>=20 && i<30)
            n_cum += nld_calc[i]*ex_bin2;
        if(i>=30 && i<40)
            n_cum += nld_calc[i]*ex_bin3;
        if(i>=40 && i<42)
            n_cum += nld_calc[i]*ex_bin4;
        if(i>=42 && i<44)
            n_cum += nld_calc[i]*ex_bin5;
        if(i>=43 && i<55)
            n_cum += nld_calc[i]*ex_bin6;
        
        /* Write Ex, ldmodel 4 temperature, n_cum, nld_calc, nld_calc (supposed to be nld_tot, with spin degeneracy) to file */
        fprintf(output89Yfile,"%7.2f %6.3f %9.2E %8.2E %8.2E ",e_calc[i],temp_calc[i],n_cum,nld_calc[i],nld_calc[i]);
        
        /* Dummy variables for next loop over spins */
        levdens = nld_calc[i];
        Ex  = e_calc[i];
        for(int j=0;j<stop_spin;j++){
            I = (double) j + start_spin;
            x = levdens*spin_distribution(Ex,I);
            fprintf(output89Yfile,"%8.2E ",x);
        }
        fprintf(output89Yfile,"\n");
        
	}
	
	
	/* Close files */	
	inputEx89Yfile.close();
	fclose(output89Yfile);

}   // END OF MAIN()

/* Constant-temperature function */
double ct_model(double Ex){
    double rho_ct = 0.;
    rho_ct = (1./T_ct)*exp((Ex - E0)/T_ct);
    cout << Ex << " " << (Ex-E0) << " " << rho_ct << endl;
    return rho_ct;
}

/* Spin function, von Egidy & Bucurescu PRC 2005/2006 */
/*double spin_distribution(double Ex, double I){
    double g_Ex = 0.;
    double U = Ex - E1;
    if (U<0.) U=0.;
    double sigma2 = 0.;

    sigma2 = 0.0146*pow(A,(5./3.))*(1. + sqrt(1. + 4.*a_par*U))/(2.*a_par);
    g_Ex = (2.*I+1.)*exp(-pow((I+0.5),2.)/(2.*sigma2))/(2.*sigma2);
    return g_Ex;
}*/

/* Spin function, von Egidy & Bucurescu PRC 2009 */
double spin_distribution(double Ex, double I){
    double g_Ex = 0.;
    double U = Ex - (0.5*Pa);
    if (U<=0.010) U=0.0010; // ! Looks strange, ask Magne (taken from counting.c)
    double sigma2 = 0.;
    
    sigma2 = 0.391*pow(A,0.675)*pow(U,0.312);
    if(sigma2 <1.) sigma2 = 1.;
    //cout << Ex << " " << U << " " << sqrt(sigma2) << endl;
    g_Ex = (2.*I+1.)*exp(-pow((I+0.5),2.)/(2.*sigma2))/(2.*sigma2);
    if(g_Ex>1.E-20) return g_Ex;
    else return 0.;
}





