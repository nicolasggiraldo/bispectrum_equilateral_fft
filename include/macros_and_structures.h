////////////////////////////////////////////////////
// HEADER WITH ALL DATA SRUCTURES FOR THE PROGRAM //
////////////////////////////////////////////////////

#define LENCHAR 500
#define X 0
#define Y 1
#define Z 2
#define INDEX(i,j,k) (k)+GV.NGRID*((j)+GV.NGRID*(i)) // Index for the C-Order
#define COMPLEXMAG(A,i) ( (A[i][0] * A[i][0]) + (A[i][1] * A[i][1]) )
#define VECTORMAG(x,y,z) sqrt( ((x)*(x)) + ((y)*(y)) + ((z)*(z)) )
#define POW2(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(X))
#define POW4(x) ((x)*(x)*(X)*(x))

/* Global variales */
struct globalVariables{
  int      NGRID;          // Number of cell in each axis
  long     NGRID2;         // Axis cell number to the square NGRID^2
  long     NGRID3;         // Total number of cells (NGRID3 = NGRID^3)
  int      GADGET_VERSION; // GADGET version of the snapshot
  double   L;              // Lenght of the simulation in Mpc.
  double   SIM_VOL;        // Volume of the simulation
  long int NP_TOT;         // Total number of particles in the simulation
  double   TOTAL_MASS;     // Total mass of all particles in the simulation
  double   RHO_MEAN;       // Mean density of ALL the simulation
  double   VOL_CELL;       // Volume of each cell
  double   H;              // Size of the cell
  double   S_KF;           // Constant to define the value of DELTA_K = S_KF * KF
  double   DELTA_K;        // Delta k for the binning
  double   KF;             // Fundamental frequency kF
  double   SHOT_NOISE;     // shotNoise of the Power Spectrum
  double   KN;             // Nyquist frequency
  char     *FILE_NAME;     // Path of the GADGET binary
  char     *OUTPUT;        // Name of the outputfile
  char     SCHEME[4];      // Scheme used for grid assignation

  /* COSMOLOGICAL PARAMETERS OF THE SIMULATION */
  double OMEGA_M0;         //Omega matter at present time
  double OMEGA_L0;         //Omega Lambda at present time
  double ZRS;              //Redshift of the simulation
  double HUBBLEPARAM;      //Hubble parameter of the simulation

  // Parameters k1,k2 and k3 assciated with sides onf the triangle
  double K1;
  //double K2;
  //double K3;

  // Power spectrum for each side and its associated errors
  double Pk1;          // Power spectrum at k1
  double Pk1_Error;    // Power spectrum error at k1
  //double Pk2;          // Power spectrum at k2
  //double Pk2_Error;    // Power spectrum error at k2
  //double Pk3;          // Power spectrum at k3
  //double Pk3_Error;    // Power spectrum error at k3

  long int Nk1;        // Number of elements in array q1
  //long int Nk2;        // Number of elements in array q2
  //long int Nk3;        // Number of elements in array q3


  double Ntri;         /* Number of counted triangles for bispectrum 
			  estimation */
  double Bk;           /* Bispectrum value to measure from the simulation 
			  snapshot in general the bispectrum is a complex                                   
			  quantity but as the avarage of density constrast                           
			  is taken the imaginary part goes to zero (see                       
			  Gil-Marin et al. 2012, J. Cosm. Astrop. Phys) */
  double Bk_Error;     // Bispectrum error at k1, k2, k3
  double Qk;           /* Reduced bispectrum, this value will be evaluated
			  with the real part of the bispectrum only */
  double Qk_Error;     // Reduced bispectrum error
  double Bk_shotnoise; // Bispectrum shotnoise
  int symmetry;        /*Symmetric factor according to the configuration 
			  if equilateral symmetry is 6 and isosceles and
			  general triangles are 2 and 1 respectively */

  /*
    
    double sumDelta2_k2; // Sum delta2 in the k1 range
    long int Nk2;        // Number of elements in array q1
    double Pk2;          // Power spectrum at k1
    double Pk2_Error;    // Power spectrum error at k1  
    double sumDelta2_k3; // Sum delta2 in the k1 range
    long int Nk3;        // Number of elements in array q1
    double Pk3;          // Power spectrum at k1
    double Pk3_Error;    // Power spectrum error at k1
  */

};
