#include "pluto.h"
#include "mutemp.h"

#define frac_Z   1.e-3   /* = N(Z) / N(H), fractional number density of 
                              metals (Z) with respect to hydrogen (H) */ 
#define frac_He  0.082   /* = N(Z) / N(H), fractional number density of 
                              helium (He) with respect to hydrogen (H) */ 

#define TRC_Z TRC+1

//double *g_Tmu_tab, *g_mu_tab;
//int g_ntabmu = 0;
//int g_oof = 0;

void LoadMuTable()
{
  FILE *fmu;
  double dummy1, dummy2;

  if (g_mu_tab == NULL) //Read mu table
  {
    print1 (" > Reading mu table from disk...\n");
    fmu = fopen("mutable.dat","r");
    if (fmu == NULL){
      print1 ("! mutemp: mutable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    g_mu_tab = ARRAY_1D(20000, double);
    g_Tmu_tab = ARRAY_1D(20000, double);

    g_ntabmu = 0;
    while (fscanf(fmu, "%lf %lf %lf %lf\n", g_Tmu_tab + g_ntabmu, &dummy1, &dummy2,
                  g_mu_tab + g_ntabmu)!=EOF) {
      g_ntabmu++;
    }
  }
}

double GetMuFromTable(double T)
{
  int    klo, khi, kmid;
  double  Tmid, dT;

  if (T > 2.e5) {return (double)g_mu_tab[g_ntabmu-1];} //Above this temperature mu is approximately constant

  klo = 0;
  khi = g_ntabmu - 1;

  if (T > g_Tmu_tab[khi])
  {
    g_oof++;//print (" ! T out of range   %12.6e\n",T);
    //    QUIT_PLUTO(1);
    return g_mu_tab[khi];
  }
  else if (T < g_Tmu_tab[klo])
  {
    g_oof++;//print (" ! T out of range   %12.6e\n",T);
    //    QUIT_PLUTO(1);
    return g_mu_tab[klo];
  }
  else
  {
    /* ----------------------------------------------
              Table lookup by binary search
       ---------------------------------------------- */
    while (klo != (khi - 1))
    {
      kmid = (klo + khi)/2;
      Tmid = g_Tmu_tab[kmid];
      if (T <= Tmid){
        khi = kmid;
      }else if (T > Tmid){
        klo = kmid;
      }
    }
    dT = g_Tmu_tab[khi] - g_Tmu_tab[klo];
    return g_mu_tab[klo]*(g_Tmu_tab[khi] - T)/dT + g_mu_tab[khi]*(T - g_Tmu_tab[klo])/dT;
  }
}

double fT(double T, void *par)
{
  struct temperature_params *p = (struct temperature_params *) par;
  p->mu = GetMuFromTable(T);
  //if (T < 2e5) {print("T: %f mu: %f\n",T,*mu);}
  //  print("T, mu, T_prs, T_mu: %f %f %f %f\n",T,*mu, prs*KELVIN/rho, T/(*mu));
  //ncall++;
  //if (ncall % 10000 == 0) {print("Process %d: ncall: %d\n", prank, ncall);}
  return p->prs*KELVIN/p->rho*p->mu - T;
}


/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
  g_minCoolingTemp = 1.e4;
  int    klo, khi, kmid;
  static int ntab;
  double  mu, T, Tmid, scrh, dT, prs;
  static double *L_tab_z0, *L_tab_z05, *L_tab_z1, *T_tab, E_cost;
  //static int niter=0;  
  FILE *fcool;
  int OutOfBounds_low = 0;
  int OutOfBounds_hi = 0;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

//Modified to read multiple tables with different metallcities
//IMPORTANT: All tables must have the exact same number of lines
//and the same temperatures

  if (T_tab == NULL){
    double dummy = 0.0;
    print1 (" > Reading tables from disk...\n");
    fcool = fopen("cooltable-z0.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-z0.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab_z0 = ARRAY_1D(7000, double);
    T_tab = ARRAY_1D(7000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                       L_tab_z0 + ntab)!=EOF) {
      ntab++;
    }
    fclose(fcool);

    fcool = fopen("cooltable-z05.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-z05.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab_z05 = ARRAY_1D(7000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", &dummy, 
                                       L_tab_z05 + ntab)!=EOF) {
      ntab++;
    }
    fclose(fcool);

    fcool = fopen("cooltable-z1.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-z1.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab_z1 = ARRAY_1D(7000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", &dummy, 
                                       L_tab_z1 + ntab)!=EOF) {
      ntab++;
    }
    fclose(fcool);
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
  }


/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

  prs = v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }



  double a, b, c;
  int count, maxiter;
  double tolerance;
  count = 0;
  tolerance = 0.05;

  double rho = v[RHO];
  //Find temperature and mu as root of mu(T)*P/(k_B*rho) - T=0 using the Brent method
  a = prs*((double)g_mu_tab[0]+tolerance)*KELVIN/rho;
  b = prs*((double)g_mu_tab[g_ntabmu-1]-tolerance)*KELVIN/rho;

  struct temperature_params par;
  par.mu = -1;
  par.rho = rho;
  par.prs = prs;
  int status = Brent(fT, &par, a, b, -1, 1.e-12, &T);
  if (status != 0)
  {
    print("! radiat: Failed to find root of temperature function on proc %d\nStatus: %d\n", prank, status);
    QUIT_PLUTO(1);
  }
  mu = par.mu;

  if (T != T){
    printf (" ! Nan found in radiat \n");
    printf (" ! rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    QUIT_PLUTO(1);
  }

  if (T < g_minCoolingTemp) { 
    rhs[RHOE] = 0.0;
    return;
  }

  int Zstatus = 0;
  double *L_tab_loZ = NULL;
  double *L_tab_hiZ = NULL;
  double Zhi = 0.0;
  double Zlo = 0.0;
  if (v[TRC_Z] > 1.0) Zstatus = 2; //Metallicity greater than max table
  else if (v[TRC_Z] > 0.3162)
  {
    Zhi = 1.0;
    Zlo = 0.3162;
    L_tab_loZ = L_tab_z05;
    L_tab_hiZ = L_tab_z0;
  }
  else if (v[TRC_Z] > 0.1)
  {
    Zhi = 0.3162;
    Zlo = 0.1;
    L_tab_loZ = L_tab_z1;
    L_tab_hiZ = L_tab_z05;
  }
  else Zstatus = 1; //Metallicity lower than min table
  double dZ = Zhi - Zlo;

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T > T_tab[khi]){
    OutOfBounds_hi++;
    //    print (" ! T out of range   %12.6e\n",T);
//    QUIT_PLUTO(1);
    if (Zstatus == 1) scrh = L_tab_z1[khi];
    else if (Zstatus == 2) scrh = L_tab_z0[khi];
    else scrh = L_tab_loZ[khi]*(Zhi - v[TRC_Z])/dZ + L_tab_hiZ[khi]*(v[TRC_Z] - Zlo)/dZ;

  }
  else if (T < T_tab[klo]){
    OutOfBounds_low++;
    //    print (" ! T out of range   %12.6e\n",T);
//    QUIT_PLUTO(1);
    if (Zstatus == 1) scrh = L_tab_z1[klo];
    else if (Zstatus == 2) scrh = L_tab_z0[klo];
    else scrh = L_tab_loZ[klo]*(Zhi - v[TRC_Z])/dZ + L_tab_hiZ[klo]*(v[TRC_Z] - Zlo)/dZ;
  }
  else
  {
    while (klo != (khi - 1))
    {
      kmid = (klo + khi)/2;
      Tmid = T_tab[kmid];
      if (T <= Tmid)
      {
        khi = kmid;
      }
      else if (T > Tmid)
      {
        klo = kmid;
      }
    }
    dT       = T_tab[khi] - T_tab[klo];
    if (Zstatus == 1) scrh = L_tab_z1[klo]*(T_tab[khi] - T)/dT + L_tab_z1[khi]*(T - T_tab[klo])/dT;
    else if (Zstatus == 2) scrh = L_tab_z0[klo]*(T_tab[khi] - T)/dT + L_tab_z0[khi]*(T - T_tab[klo])/dT;
    else
    {
      double L_loZ = L_tab_loZ[klo]*(T_tab[khi] - T)/dT + L_tab_loZ[khi]*(T - T_tab[klo])/dT;
      double L_hiZ = L_tab_hiZ[klo]*(T_tab[khi] - T)/dT + L_tab_hiZ[khi]*(T - T_tab[klo])/dT;
      scrh     = L_loZ*(Zhi - v[TRC_Z])/dZ + L_hiZ*(v[TRC_Z] - Zlo)/dZ;
    }
  }

  rhs[RHOE] = -scrh*v[RHO]*v[RHO];
  rhs[RHOE] *= E_cost*UNIT_DENSITY*UNIT_DENSITY/(mu*mu*CONST_mp*CONST_mp);

//  print("T, Z, L, rhs: %6.2e %f %6.2e %6.2e\n", T, v[TRC_Z], scrh, rhs[RHOE]);


#ifdef CELLINFO_VERBOSE
  if (OutOfBounds_low > 0)  print("! radiat: %d cells on proc %d had temperatures below the minimum in table\n", OutOfBounds_low, prank);
  if (OutOfBounds_hi > 0)  print("! radiat: %d cells on proc %d had temperatures above the maximum in table\n", OutOfBounds_hi, prank);
#endif
}
#undef T_MIN
/* ******************************************************************* */
double MeanMolecularWeight (double *V)
/*
 *
 *
 *
 ********************************************************************* */
{
  //This assumes fully ionised gas, to properly calculate the mean molecular weight for a given temperature call GetMuFromTable instead
  return (0.59);//v[TRC_MU]);
/*
  return  ( (A_H + frac_He*A_He + frac_Z*A_Z) /
            (2.0 + frac_He + 2.0*frac_Z - 0.0));
*/
}



