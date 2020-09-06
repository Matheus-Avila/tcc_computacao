#include "IS_Model.h"

/******************************************************************************
 * 
 * IS_Model - Imunne System Model 
 * 
 * This software simmulates the behavior of the main cells of the immune
 * system both innate and acquired.
 * 
 * Based on the works from Pigozzo (2011) and Marchuk (1997).
 * 
 * Date modified : 06/12/2013.
 * 
 * Recquires: 'IS_Model.h'.
 * 
 * Outputs : 
 *  
 *          '*.csv' files for time steps chosen for each PDE in the model.
 *          'L.dat' containing the averages of cells in the tissue.
 *          'T.dat', 'B.dat', 'P.dat' containing these cells concentrations 
 *                                    over time.
 * 
 * Use-me : 
 *  
 *          1. Instantiate an object of IS_Model class:
 * 
 *             IS_Model* model = new IS_Model();
 * 
 *          2. Call solve() method:
 *  
 *             model->solve();
 * 
 * 
 * Obs: Considers logistic growth of bacteria.
 * 
 *  
 ******************************************************************************/

using namespace std;

/**
 * Greetings from the simulator
 */
std::string IS_Model::Header()
{
    int size = 40;
    std::string returnstring;
    std::ostringstream sstream;
    sstream << a0;
    std::string a0str = sstream.str();

    for(int i = 0; i < size; i++) returnstring += "*";
    returnstring += "\n* Begin of simulation SI3D. \n";
    returnstring += "*\n* Initial bacteria = "+a0str+".\n*\n";
    for(int i = 0; i < size; i++) returnstring += "*";
    returnstring += "\n";
    return returnstring;
}

/**
 * Goodbye messagem from the simulator
 */
std::string IS_Model::Footer(long int t){
    std::string returnstring;
    int size = 40;
    std::ostringstream sstream1, sstream2;
    sstream1 << t;
    sstream2 << days;
    std::string tstr = sstream1.str();
    std::string daysstr = sstream2.str();

    returnstring += "The end. \n ...of simulation!";
    returnstring += "\n"+tstr+" time steps for "+daysstr+" days.\n";
    for(int i = 0; i < size; i++) returnstring += "*";
    sstream1 << A_T;
    std::string atstr = sstream1.str();

    returnstring += "\nBacteria in the end : "+atstr+" \n";
    for(int i = 0; i < size; i++) returnstring += "*";
    returnstring += "\n";
    return returnstring;
}

/**
 * Setters
 */
void IS_Model::setSaveFiles(int sf){
    this->saveFiles = sf;
}
void IS_Model::setSimulationCase(int sc){
    this->simCase = sc;
}

/**
* Constructor
*/
IS_Model::IS_Model(){

}

/**
* Set conditions and parameter values
*/
void IS_Model::initialize(){

  /** 
   * 0 - simulates coupled model, 
   * 1 - simulates only antigen diffusion,
   * 2 - simulates only innate response,
   * 3 - complete model without diffusion.
   */  
  simCase   = 0;  
  /** 
   * 0 - saves only edos files, 
   * 1 - saves all files.
   */
  saveFiles = 1;
  /**
   * number of days simulated
   */
  days      =  1; //change the time interval here
  /**
   * number of files saved
   */
  points    = 10;
  /**
   * output directory
   */
  dir       = (char *) "output/";
  /**
   * 0 - contact with lymph vessels only on one border, 
   * 1 - homogeneous contact with lymph vessels.
   * 2 - contact with lymph vessels given by function.
   */
  lnv       = 2;
  /**
   * 0 - contact with blood vessels only on one border,
   * 1 - homogeneous contact with blood vessels,  
   * 2 - contact with blood vessels given by function.
   */
  bv        = 2;
  /**
   * each (5/(pow(10,6)) = 0.0000002 days or 0,01728 secs
   */
  deltaT     = pow(10,-4); //day
  /**
   * each 1000000 iterations represents 1 day
   */
  iterPerDay = pow(10,4);
  /**
   * each deltaX represents 100 micrometers ((1 × 10^-6 m)), a cell
   * has 1000 cubic micrometers, each discretized space has 1.000.000 cubic micrometers,
   * leading to 1000 cells for each space
   */
  deltaX     = 0.1; //mm
  deltaY     = 0.1; //mm
  deltaZ     = 0.1; //mm
  
  tol        = pow(10,-6);
  //  Inicializado no main e passado pelo solve
  m0         = 0.0;
  a0         = 2.0;
  th0        = 0.0;
  b0         = 0.0;
  p0         = 0.0;

  //immune system initial condition csnnot be 0!!!!!!!!!!!!!!!!!!!
  TNF_0       =  0.04;
  IL_10_0      = 0.00027;
  IL_6_0      = 0.04;
  IL_8_0      = 0.02;


  //CREA regulation constants
  n_610 = 34.8/1000; //pg/mm^3
  n_66 = 560/1000; // pg/mm^3
  n_6TNF = 185/1000; // pg/mm^3
  n_810 = 17.4/1000; // pg/mm63
  n_8TNF = 185/1000; // pg/mm^3
  n_106 = 560/1000; // pg/mm^3
  n_TNF10 = 17.4/1000; // pg/mm^3
  n_TNF6 = 560/1000; // pg/mm^3
  n_M10 = 4.35/1000; // pg/mm^3
  n_MTNF = 0.1; // pg/mm^3

  h_106 = 3.68;
  h_6TNF = 2;
  h_66 = 1;
  h_610 = 4;
  h_8TNF = 3;
  h_810 = 1.5;
  h_TNF10 = 3;
  h_TNF6 = 2;
  h_M10 = 0.3;
  h_MTNF = 3.16;


  t_estrela  = 8.4 * pow(10,-3);// * MOL;
  b_estrela  = 8.4 * pow(10,-4);// * MOL;
  p_estrela  = 8.4 * pow(10,-6);// * MOL;
  m_estrela  = 4.0;// * MOL
  f0         = 0;// * MOL;
  f_estrela  = 0;// * MOL;

  MA_T    = 1.0;
  MA_L    = 0.0;
  MR_T    = m_estrela;
  Th      = th0;
  B       = b0;
  P       = p0;
  F_T     = f0;
  F_L     = f0;//f_estrela;
  A_T     = a0/1000;
  TNF_alpha_T = TNF_0;
  IL_6_T = IL_6_0;
  IL_8_T = IL_8_0;
  IL_10_T = IL_10_0;
  d_a        = 0.00037;       // antigen diffusion (Haessler) mm^3/day
  d_mr       = 0.0432;        //resting macrophage diffusion (estimated) mm^3/day
  d_ma       = 0.3;           //macrophage diffusion (estimated) mm^3/day
  d_f        = 0.016;          //antibody diffusion (estimated) mm^3/day


//TNF
  k_TNFM = 0.6 * 24/(1000); //pg/mm^3*day* # of cells
  k_TNF = 24; // day^-1
  q_TNF = 4; //pg/mm^3
//IL_6
  k_6m = 0.81 * 24/1000; //pg/mm^3*day* # of cells
  k_6TNF = 0.81 * 24 /1000; //pg/mm^3*day* # of cells
  k_6 = 4.64*24;//day^-1
  q_IL6 = 4; //pg/mm^3
//IL_8
  k_8m = 0.56 *24/1000;//pg/mm^3*day* # of cells
  k_8TNF = 0.56*24/1000; //pg/mm^3*day* # of cells
  k_8 = 4.64*24; //day^-1
  q_IL8 = 2; //pg/mm^3
//IL_10
  k_10m = 0.0191*24/1000; //pg/mm^3*day* # of cells
  k_106 = 0.0191 *24/1000; //pg/mm^3*day* # of cells
  k_10 = 0.5*24; //day^-1
  q_IL10 = 0.27; // pg/mm^3


  beta_A     = 2.0;           //bacteria replication 1/day
  k_A        = 50.0;          //maximum capacity cell/mm^3
  m_A        = 0.1;            //bacteria natural decay 1/day
  m_Mr       = 0.033;          //resting macrophage natural decay 1/day
  m_Ma       = 0.07;           //activated macrophage natural decay 1/day

  k_m = 0.0414 *24; // 1/day
  k_MTNF = 8.65 * 24; // 1/day

  
  gamma_ma   = 8.30*pow(10,-2); //macrophage activation mm^3/(cell * day)
  
  lambda_ma  = 5.98*pow(10,-2); //activated macrophage fagocitosis rate mm^3/(cell * day)
  lambda_mr  = 5.98*pow(10,-3); //resting macrophage fagocitosis rate mm^3/(cell * day)
  lambda_afma= 7.14*pow(10,-2); //activated macrophage fagocitosis rate for opsonized antigen mm^3/(cell * day)
  lambda_afmr= 1.66*pow(10,-3); //resting macrophage fagocitosis rate for opsonized antigen mm^3/(cell * day)

  b_th       = 1.7*pow(10,-2); //th2 stimuli mm^3/(cell*day)
  b_p        = 1.*pow(10,5);  //th2 expenditure to stimulate b mm^6/(cell*day)
  b_pb       = 6.02*pow(10,3);//b stimuli mm^6/(cell*day)
  b_pp       = 2.3*pow(10,6); //b stimuli while describes plasma cell mm^3/(cell*day)

  ro_t       = 2.0;            //th2 descendents -
  ro_b       = 16.0;           //b descendents Cell/mm^3
  ro_p       = 3.0;            //p descendents -
  ro_f       = 5.1*pow(10,4); //antibody release -

  alpha_t    = 0.01;           //th2 natural decay 1/day
  alpha_b    = 1.0;            //b natural decay 1/day
  alpha_p    = 5.0;            //plasma decay 1/day
  alpha_f    = 0.43;          //antibody migration 1/day
  alpha_Ma   = 0.001;// * pow(10,1);  //migration (estimated) 1/day
  alpha_mr   = 4.0;            //resting macrophage source coefficient 1/day
  
  /**
   * Initial Conditions
   */
  /*for(int x = 0; x < Xspace; x++) 
  {
    for(int y = 0; y < Yspace; y++) 
    {
      for(int z = 0; z < Zspace; z++) 
      {
        if (simCase == 3)
        {
          A[0][x][y][z] = a0/10000;
        }
        else
        {
          //bacteria only in the center of the cubic domain
          if ((x > (0.2*Xspace)&&( x < (0.7*Xspace)))
            && (y > (0.2*Yspace)&&( y < (0.7*Yspace)))
            && (z > (0.2*Zspace)&&( z < (0.7*Zspace)))) 
          {
            A[0][x][y][z] = a0;///IC_SPACE;
          } 
          else
          {
            A[0][x][y][z] = 0.0;
          }   
        }
        MR[0][x][y][z]  = m_estrela;
        MA[0][x][y][z]  = 1.0;
        F[0][x][y][z]   = f0;//SPACE;
        TNF_alpha[0][x][y][z] = TNF_0;
        
        IL_6[0][x][y][z] = IL_6_0;
        IL_8[0][x][y][z] =IL_8_0;
        IL_10[0][x][y][z] =IL_10_0;
        
      }
    }
  }*/
}

/**
 * Updates the current results to position 1
 * and sets the previous results to zero
 */
void IS_Model::update(double vec[][Xspace][Yspace][Zspace])
{
  for(int x = 0; x < Xspace; x++) 
  {
    for(int y = 0; y < Yspace; y++) 
    {
      for(int z = 0; z < Zspace; z++)
      {
        vec[0][x][y][z] = vec[1][x][y][z];
        /**
         * In case it is necessary to keep the previous
   * just comment the line below
   */
  //vec[1][x][y][z] = 0.;
      }
    }
  }
}

/**
 * Calculates the laplacian for given value and position
 */
double IS_Model::laplacian(double vec[][Xspace][Yspace][Zspace], int x, int y, int z){
  double resX = 0, resY = 0, resZ = 0;
  // same boundary condition to every equation
  if(x == 0) 
  {
    resX = (vec[0][x+1][y][z] - vec[0][x][y][z])/(deltaX*deltaX);
  } 
  else if(x == Xspace-1 ) 
  {
    resX = (vec[0][x-1][y][z] - vec[0][x][y][z])/(deltaX*deltaX);
  } 
  else 
  {//dentro do dominio mas fora da extremidade
    resX = (vec[0][x+1][y][z] -2 * vec[0][x][y][z] + vec[0][x-1][y][z])/(deltaX*deltaX);
  } // X
  if(y == 0) 
  {
    resY = (vec[0][x][y+1][z] - vec[0][x][y][z])/(deltaY*deltaY);
  } 
  else if( y == Yspace-1) 
  {
    resY = (vec[0][x][y-1][z] - vec[0][x][y][z])/(deltaY*deltaY);
  } 
  else 
  {
    resY = (vec[0][x][y+1][z] -2 * vec[0][x][y][z] + vec[0][x][y-1][z])/(deltaY*deltaY);
  } //Y
  if(z == 0) 
  {
    resZ = (vec[0][x][y][z+1] - vec[0][x][y][z])/(deltaZ*deltaZ);
  } 
  else if( z == Zspace-1) 
  {
    resZ = (vec[0][x][y][z-1] - vec[0][x][y][z])/(deltaZ*deltaZ);
  } 
  else 
  {
    resZ = (vec[0][x][y][z+1] -2 * vec[0][x][y][z] + vec[0][x][y][z-1])/(deltaZ*deltaZ);
  }//Z
  return resX+resY+resZ;
}

/**
 * Tests if the file could be created and exits the simulation if
 * there was an error
 */
int IS_Model::checkFile(FILE* theFile)
{
  if (theFile==NULL)
  {
    cout << "Error opening file!!!\n Make sure the path is correct! \n";
    return 1;
  }
  return 0;
}

/**
 * Calculates integrals of cells in the tissue and return the value as
 * a pointer 
 */
int IS_Model::calcIntegral(double vec[][Xspace][Yspace][Zspace], double *V){

  for(int x = 0; x < Xspace; x++) 
  {
    for(int y = 0; y < Yspace; y++) 
    {
      for(int z = 0; z < Zspace; z++) 
      {
        if (vec[0][x][y][z]>0.0) *V += vec[0][x][y][z];
      }
    }
  }
  if (*V > 0.0) *V = (*V/(SPACE)); 
  else *V = 0.0; //V
  return 0;
}

/**
 * For activated macrophages the integral is calculated considering only 
 * the cells in contact with lymph vessels
 */
int IS_Model::calcIntegral_lv(double vec[][Xspace][Yspace][Zspace], double *V){

  for(int x = 0; x < Xspace; x++) 
  {
    for(int y = 0; y < Yspace; y++) 
    {
      for(int z = 0; z < Zspace; z++) 
      {
        if (((lnv==0)&&(x==0))||(lnv==1)||((lnv==2)&&(is_lnvase(x,y,z))))
        {
          if (vec[0][x][y][z]>0.0) *V += vec[0][x][y][z];
        }
      }
    }
  }
  if (*V > 0.0) *V = (*V/(SPACE)); else *V = 0.0;
  return 0;
}

/**
 * for antibodies consider only cells in contact with blood vessels
 */
int IS_Model::calcIntegral_bv(double vec[][Xspace][Yspace][Zspace], double *V){

  for(int x = 0; x < Xspace; x++)
  {
    for(int y = 0; y < Yspace; y++) 
    {
      for(int z = 0; z < Zspace; z++) 
      {
        if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z))))
        {
          if (vec[0][x][y][z]>0.0) *V += vec[0][x][y][z];
        }
      }
    }
  }
  if (*V > 0.0) *V = (*V/(SPACE)); 
  else *V = 0.0;
  return 0;
}

/**
 * returns 1 if the point is a blood vase and 0 if it is not
 */
int IS_Model::is_bvase(int x, int y, int z)//do not understand
{
  //if(((x >= 0)&&(x <= 1))||((x>=4)&&(x<=5))||((x>=8)&&(x<=9)))
  //  if(((z >= 0)&&(z <= 1))||((z>=4)&&(z<=5))||((z>=8)&&(z<=9)))
 if(((x >= 0)&&(x <= 1))||((x>=8)&&(x<=9)))
    if(((z >= 0)&&(z <= 1))||((z>=8)&&(z<=9)))  //no error here?
      return 1;
  
  return 0;     
}

/**
 * returns 1 if the point is a lymph vase and 0 if it is not
 */
int IS_Model::is_lnvase(int x, int y, int z)
{
  if(((x >= 2)&&(x <= 3))||((x>=6)&&(x<=7)))
    if(((z >= 0)&&(z <= 1))||((z>=4)&&(z<=5)))
      return 1;
  
  return 0;     
}


/******************************************************************************
* Solve model equations
*******************************************************************************/
double* IS_Model::solve(double* parameters)
{
  int i       = 0;
  long int t  = 0;

  //set initial conditions
  //Ainda precisa remover as inicializações do initialize()
  initialize();

  a0 = parameters[0];  
  
  TNF_0 =  parameters[1];
  
  IL_10_0 = parameters[2];
  
  IL_6_0 = parameters[3];
  
  IL_8_0 = parameters[4];


  n_610 = parameters[5]; //pg/mm^3
  
  n_66 = parameters[6]; // pg/mm^3
  
  n_6TNF = parameters[7]; // pg/mm^3
  
  n_810 = parameters[8]; // pg/mm63
  
  n_8TNF = parameters[9]; // pg/mm^3
  
  n_106 = parameters[10]; // pg/mm^3
  
  n_TNF10 = parameters[11]; // pg/mm^3
  
  n_TNF6 = parameters[12]; // pg/mm^3
  
  n_M10 = parameters[13]; // pg/mm^3
  
  n_MTNF = parameters[14]; // pg/mm^3
  
  h_106 = parameters[15];
  
  h_6TNF = parameters[16];
  
  h_66 = parameters[17];
  
  h_610 = parameters[18];
  
  h_8TNF = parameters[19];
  
  h_810 = parameters[20];
  
  h_TNF10 = parameters[21];
  
  h_TNF6 = parameters[22];
  
  h_M10 = parameters[23];
  
  h_MTNF = parameters[24];
  
  
  A_T     = parameters[25];
  

  TNF_alpha_T = parameters[26];
  
  IL_6_T = parameters[27];
  
  IL_8_T = parameters[28];
  
  IL_10_T = parameters[29];
  
  d_a        = parameters[30];       // antigen diffusion (Haessler) mm^3/day
  
  d_mr       = parameters[31];        //resting macrophage diffusion (estimated) mm^3/day
  
  d_ma       = parameters[32];           //macrophage diffusion (estimated) mm^3/day
  
  d_f        = parameters[33];          //antibody diffusion (estimated) mm^3/day
  


//TNF
  k_TNFM = parameters[34]; //pg/mm^3*day* # of cells
  
  k_TNF = parameters[35]; // day^-1

  q_TNF = parameters[36]; //pg/mm^3
//IL_6
  k_6m = parameters[37]; //pg/mm^3*day* # of cells
  
  k_6TNF = parameters[38]; //pg/mm^3*day* # of cells
  
  k_6 = parameters[39];//day^-1
  
  q_IL6 = parameters[40]; //pg/mm^3
    
//IL_8
  k_8m = parameters[41];//pg/mm^3*day* # of cells
  
  k_8TNF = parameters[42]; //pg/mm^3*day* # of cells
  
  k_8 = parameters[43]; //day^-1
  
  q_IL8 = parameters[44]; //pg/mm^3
    
//IL_10
  k_10m = parameters[45]; //pg/mm^3*day* # of cells
  
  k_106 = parameters[46]; //pg/mm^3*day* # of cells
  
  k_10 = parameters[47]; //day^-1
  
  q_IL10 = parameters[48]; // pg/mm^3
  

  k_m = parameters[49]; // 1/day
  
  k_MTNF = parameters[50]; // 1/day

for(int x = 0; x < Xspace; x++) 
  {
    for(int y = 0; y < Yspace; y++) 
    {
      for(int z = 0; z < Zspace; z++) 
      {
        if (simCase == 3)
        {
          A[0][x][y][z] = a0/10000;
        }
        else
        {
          //bacteria only in the center of the cubic domain
          if ((x > (0.2*Xspace)&&( x < (0.7*Xspace)))
            && (y > (0.2*Yspace)&&( y < (0.7*Yspace)))
            && (z > (0.2*Zspace)&&( z < (0.7*Zspace)))) 
          {
            A[0][x][y][z] = a0;///IC_SPACE;
          } 
          else
          {
            A[0][x][y][z] = 0.0;
          }   
        }
        MR[0][x][y][z]  = m_estrela;
        MA[0][x][y][z]  = 1.0;
        F[0][x][y][z]   = f0;//SPACE;
        TNF_alpha[0][x][y][z] = TNF_0;
        
        IL_6[0][x][y][z] = IL_6_0;
        IL_8[0][x][y][z] =IL_8_0;
        IL_10[0][x][y][z] =IL_10_0;
        
      }
    }
  }

  //print program header
  cout << Header();

  //alloc memory for filename
  //char *fileName = (char *)malloc(20*sizeof(char));
  
   // sprintf(fileName, "%s%s", dir, "L.csv");
   // datamatlabL = fopen(fileName, "w");

    //CREA made
   // sprintf(fileName, "%s%s", dir, "Cytokines.csv");
   // datamatlabCytokine = fopen(fileName, "w");
    
    //check valid dir
    //if (checkFile(datamatlabL)) return 1;
   // if (checkFile(datamatlabCytokine)) return 1;
 

    //sprintf(fileName, "%s%s", dir, "T.dat");
    //datamatlabT = fopen(fileName, "w");

    //sprintf(fileName, "%s%s", dir, "B.dat");
    //datamatlabB = fopen(fileName, "w");

    //sprintf(fileName, "%s%s", dir, "P.dat");
    //datamatlabP = fopen(fileName, "w");

/**
 * begin time loop
 */
  do
  {

    if (t==0) cout << "Calculating...\n";

    //int value = ((int)iterPerDay*days)/points;

    /*if(t%value == 0) 
    {
      cout << "Saving files : iteration ..."<< t << "\n";

      if (saveFiles)
      {//check before saving all files
        sprintf(fileName,"%sA_%ld.csv", dir,t);
        datamatlabA = fopen(fileName, "w");
        sprintf(fileName,"%sMr_%ld.csv",dir,t);
        datamatlabMr = fopen(fileName, "w");
        sprintf(fileName,"%sMa_%ld.csv",dir,t);
        datamatlabMa = fopen(fileName, "w");
        sprintf(fileName,"%sF_%ld.csv",dir,t);
         datamatlabF = fopen(fileName, "w");

         //add sth

          sprintf(fileName,"%sTNF_alpha_%ld.csv", dir,t);
          datamatlabTNF_alpha = fopen(fileName, "w");
          
          //sprintf(fileName,"%sIL_8_%ld.csv", dir,t);
          //datamatlabIL_8 = fopen(fileName, "w");
          //sprintf(fileName,"%sIL_10_%ld.csv", dir,t);
          //datamatlabIL_10 = fopen(fileName, "w");
          //sprintf(fileName,"%sIL_6_%ld.csv", dir,t);
          //datamatlabIL_6 = fopen(fileName, "w"); //If u want the detailed data just unquote all related parts,
        fprintf(datamatlabT, "%ld, %.2E \n", t, Th);
        fprintf(datamatlabB, "%ld, %.2E \n", t, B);
        fprintf(datamatlabP, "%ld, %.2E \n", t, P);

        fprintf(datamatlabL, "%ld, %.2E, %.2E, %.2E, %.2E, %.2E, %.2E \n", 
                t, MA_T, F_T, MA_L, F_L, A_T, MR_T);

        fprintf(datamatlabCytokine, "%ld, %.6E, %.6E, %.6E,%.6E\n",t, TNF_alpha_T, IL_8_T, IL_10_T, IL_6_T);
        for(int x = 0; x < Xspace; x++) 
        {
          for(int y = 0; y < Yspace; y++) 
          {
            for(int z = 0; z < Zspace; z++) 
            {
              if( (x+1 == Xspace && y+1 == Yspace && z+1==Zspace) ) 
              {
                fprintf(datamatlabA, "%d, %d, %d, %E", x, y, z, A[0][x][y][z]);
                fprintf(datamatlabMr, "%d, %d, %d, %E", x, y, z, MR[0][x][y][z]);
                fprintf(datamatlabMa, "%d, %d, %d, %E", x, y, z, MA[0][x][y][z]);
                fprintf(datamatlabF, "%d, %d, %d, %E", x, y, z, F[0][x][y][z]);

                fprintf(datamatlabTNF_alpha, "%d, %d, %d, %E", x, y, z, TNF_alpha[0][x][y][z]);
                //fprintf(datamatlabIL_8, "%d, %d, %d, %E", x, y, z, IL_8[0][x][y][z]);
                //fprintf(datamatlabIL_10, "%d, %d, %d, %E", x, y, z, IL_10[0][x][y][z]);
                //fprintf(datamatlabIL_6, "%d, %d, %d, %E", x, y, z, IL_6[0][x][y][z]);
              } 
              else 
              {
                fprintf(datamatlabA,  "%d, %d, %d, %E\n", x, y, z, A[0][x][y][z]);
                fprintf(datamatlabMr, "%d, %d, %d, %E\n", x, y, z, MR[0][x][y][z]);
                fprintf(datamatlabMa, "%d, %d, %d, %E\n", x, y, z, MA[0][x][y][z]);
                 fprintf(datamatlabF, "%d, %d, %d, %E\n", x, y, z, F[0][x][y][z]);  
                fprintf(datamatlabTNF_alpha, "%d, %d, %d, %E\n", x, y, z, TNF_alpha[0][x][y][z]);
                //fprintf(datamatlabIL_8,"%d, %d, %d, %E\n", x, y, z, IL_8[0][x][y][z]);
                //fprintf(datamatlabIL_10, "%d, %d, %d, %E\n", x, y, z, IL_10[0][x][y][z]);
                //fprintf(datamatlabIL_6, "%d, %d, %d, %E\n", x, y, z, IL_6[0][x][y][z]);
              }
            }
          }
        }
        fclose(datamatlabA);
        fclose(datamatlabMr);
        fclose(datamatlabMa);
        fclose(datamatlabF);
        fclose(datamatlabTNF_alpha);
        //fclose(datamatlabIL_8);
        //fclose(datamatlabIL_10);
        //fclose(datamatlabIL_6); 
      }
      else
      {
        fprintf(datamatlabT, "%ld, %.2E \n", t, Th);
        fprintf(datamatlabB, "%ld, %.2E \n", t, B);
        fprintf(datamatlabP, "%ld, %.2E \n", t, P);
        fprintf(datamatlabL, "%ld, %.2E, %.2E, %.2E, %.2E, %.2E, %.2E\n", t, MA_T, F_T, MA_L, F_L, A_T, MR_T);  

        //CREA made
        fprintf(datamatlabCytokine, "%ld, %.2E, %.2E, %.2E, %.2E\n", t,TNF_alpha_T, IL_8_T, IL_10_T, IL_6_T);
  
//      fprintf(datamatlabL, "%ld %.2E %.2E %.2E %.2E %.2E %.2E \n", t,
                //MA[0][0][0][0], F[0][0][0][0], MA_L, F_L, A[0][0][0][0],
                //MR[0][0][0][0]);
      }
    }
  */
    //integral
    //cout << "Solve integrals. ";
    if (t > 0 && simCase!=3)
    {
      MA_T = MR_T = F_T = A_T = IL_6_T = TNF_alpha_T = IL_8_T = IL_10_T = 0.0;
      if (calcIntegral_lv(MA, &MA_T)!=0)
      {
        cout << "Something went wrong with the integral!!! \n";
        double* saida;
        return saida;
      }      
      calcIntegral(MR, &MR_T);
      calcIntegral(MA, &MA_T);
      calcIntegral_bv(F, &F_T);
      calcIntegral(A, &A_T);

      calcIntegral(TNF_alpha, &TNF_alpha_T);
      calcIntegral(IL_6, &IL_6_T);
      calcIntegral(IL_8, &IL_8_T);
      calcIntegral(IL_10, &IL_10_T);

      //cout << "T = " << t << "TNFalpha: " << TNF_alpha_T << "\n";
      //cout << "IL6: "  << IL_6_T << "\n";
      //cout << "IL8: "  << IL_8_T << "\n";
      //cout << "IL10: " <<  IL_10_T << "\n";

      //calcIntegral(CH, &CH_T);
      //calcIntegral(CA, &CA_T);
      //calcIntegral(D, &D_T);
    }

//*****************************************************************************
    //Complete model without diffusion
    if (simCase == 3)
    {// Add our equations here without diffusion.
        A[0][0][0][0] = ( beta_A*A[0][0][0][0]*(1-(A[0][0][0][0]/k_A))
          - ( lambda_mr*MR[0][0][0][0]*A[0][0][0][0])
          - (lambda_ma* MA[0][0][0][0]* A[0][0][0][0])
          - (lambda_afma*F[0][0][0][0]*A[0][0][0][0]*MA[0][0][0][0])
          - (lambda_afmr*F[0][0][0][0]*A[0][0][0][0]*MR[0][0][0][0])
          - m_A * A[0][0][0][0]) * deltaT + A[0][0][0][0];

        MR[0][0][0][0] = ((- m_Mr * MR[0][0][0][0])
           - (gamma_ma * 
            (k_m + k_MTNF * pow(TNF_alpha[0][0][0][0], h_MTNF)/(pow(n_MTNF, h_MTNF) +pow(TNF_alpha[0][0][0][0], h_MTNF)))
             * pow(n_M10, h_M10)/(pow(n_M10, h_M10) +pow(IL_10[0][0][0][0], h_M10)) * 
            MR[0][0][0][0] * A[0][0][0][0])       
           + alpha_mr * (m_estrela - MR[0][0][0][0])
                         ) * deltaT + MR[0][0][0][0];

        MA[0][0][0][0] = ((-m_Ma * MA[0][0][0][0])
           + (gamma_ma * 
            (k_m + k_MTNF * pow(TNF_alpha[0][0][0][0], h_MTNF)/(pow(n_MTNF, h_MTNF) +pow(TNF_alpha[0][0][0][0], h_MTNF)))
             * pow(n_M10, h_M10)/(pow(n_M10, h_M10) +pow(IL_10[0][0][0][0], h_M10)) * 
            MR[0][0][0][0] * A[0][0][0][0])
           - alpha_Ma * (MA_T - MA_L)
            ) * deltaT + MA[0][0][0][0];

        F[0][0][0][0] = (- (lambda_afma * F[0][0][0][0]* A[0][0][0][0]*MA[0][0][0][0])
          - (lambda_afmr*F[0][0][0][0]*A[0][0][0][0]*MR[0][0][0][0])
           - (alpha_f * (F_T - F_L))
            )* deltaT + F[0][0][0][0];

        IL_6[0][0][0][0] = ((k_6m + k_6TNF * pow(TNF_alpha[0][0][0][0], h_6TNF)/(pow(n_6TNF, h_6TNF) +pow(TNF_alpha[0][0][0][0], h_6TNF)))
           * (pow(n_66, h_66)/(pow(n_66, h_66) +pow(IL_6[0][0][0][0], h_66))) * pow(n_610, h_610)/(pow(n_610, h_610) +pow(IL_10[0][0][0][0], n_610))
            * MA[0][0][0][0]
            - k_6 * (IL_6[0][0][0][0] - MA[0][0][0][0]/ MA[0][0][0][0] * MR[0][0][0][0]/MR[0][0][0][0] * q_IL6)) * deltaT + IL_6[0][0][0][0];

        TNF_alpha[0][0][0][0] = (k_TNFM * pow(n_TNF6, h_TNF6)/(pow(n_TNF6, h_TNF6) +pow(IL_6[0][0][0][0], h_TNF6))
         * pow(n_TNF10,h_TNF10)/(pow(n_TNF10, h_TNF10) +pow(IL_10[0][0][0][0], h_TNF10)) * MA[0][0][0][0]
           - k_TNF * (TNF_alpha[0][0][0][0] - MR[0][0][0][0]/MR[0][0][0][0] * q_TNF)) * deltaT + TNF_alpha[0][0][0][0];

        IL_10[0][0][0][0] = ((k_10m + k_106 * pow(IL_6[0][0][0][0], h_106)/(pow(n_106, h_106) +pow(IL_6[0][0][0][0], h_106)))
         * MA[0][0][0][0] 
          - (k_10 * (IL_10[0][0][0][0] - MA[0][0][0][0]/ MA[0][0][0][0] * MR[0][0][0][0]/MR[0][0][0][0] * q_IL10)) 
              ) * deltaT + IL_10[0][0][0][0];

        IL_8[0][0][0][0] = ((k_8m + k_8TNF * pow(TNF_alpha[0][0][0][0], h_8TNF)/(pow(TNF_alpha[0][0][0][0], h_8TNF) +pow(n_8TNF, h_8TNF)))
         * pow(n_810, h_810)/(pow(n_810, h_810) +pow(IL_10[0][0][0][0], h_810)) * MA[0][0][0][0]
          - k_8 * (IL_8[0][0][0][0] - MA[0][0][0][0]/ MA[0][0][0][0] * MR[0][0][0][0]/MR[0][0][0][0] * q_IL8)
              ) * deltaT + IL_8[0][0][0][0]; 


        MA_L = ( alpha_Ma * (MA_T - MA_L)) * deltaT + MA_L;
        if (MA_L < 0.0) MA_L = 0.0;

        Th = (b_th*(ro_t*Th*MA_L -Th*MA_L) -b_p*MA_L*Th*B
               + alpha_t*(t_estrela - Th)) * deltaT + Th;
        if (Th < 0.0) Th = t_estrela;

        B = (b_pb*(ro_b*Th*MA_L-Th*MA_L*B)
              + alpha_b*(b_estrela - B)) * deltaT + B;
        if (B < 0.0) B = b_estrela;

        P = (b_pp*(ro_p*Th*MA_L*B) + alpha_p*(p_estrela - P)) * deltaT + P;
        if (P < 0.0) P = p_estrela;

        F_L = (ro_f*P + alpha_f*(F_T-F_L)) * deltaT + F_L;
        if (F_L < 0.0) F_L = f_estrela;
    }
    
//*****************************************************************************

    else
    {
      //Solve ODEs    
      if(simCase==0)
      {
        MA_L = ( alpha_Ma * (MA_T - MA_L)) * deltaT + MA_L;
        if (MA_L < 0.0) MA_L = 0.0;

        Th = (b_th*(ro_t*Th*MA_L -Th*MA_L) -b_p*MA_L*Th*B
               + alpha_t*(t_estrela - Th)) * deltaT + Th;
        if (Th < 0.0) Th = t_estrela;

        B = (b_pb*(ro_b*Th*MA_L-Th*MA_L*B)
              + alpha_b*(b_estrela - B)) * deltaT + B;
        if (B < 0.0) B = b_estrela;

        P = (b_pp*(ro_p*Th*MA_L*B) + alpha_p*(p_estrela - P)) * deltaT + P;
        if (P < 0.0) P = p_estrela;

        F_L = (ro_f*P + alpha_f*(F_T-F_L)) * deltaT + F_L;
        if (F_L < 0.0) F_L = f_estrela;
      }

    //Solve PDEs
        // Add our equations with diffusion here.
      for(int x = 0; x < Xspace; x++)
      {
        for(int y = 0; y < Yspace; y++) 
        {
          for(int z = 0; z < Zspace; z++) 
          {
    
//*****************************************************************************
          //Simulates only antigen diffusion
            if(simCase==1)
            {
    
      //Antigenos
              A[1][x][y][z] = ( beta_A*A[0][x][y][z]*(1-(A[0][x][y][z]/k_A)) 
                + (d_a * laplacian(A,x,y,z))
                - m_A * A[0][x][y][z]) * deltaT + A[0][x][y][z];

              if(A[1][x][y][z] != A[1][x][y][z]) 
              {
                //cout << "A\t(NaN)-> i:"<<i<<"\t-> ("<< x << y << z << ")" << A[1][x][y][z] << "\n";
      //  return 1;
              }
//*****************************************************************************
    //Simulates only innate response
            }
            else if (simCase==2)
            {

            //Antigenos
              A[1][x][y][z] = ( beta_A*A[0][x][y][z]*(1-(A[0][x][y][z]/k_A))
                - ( lambda_mr*MR[0][x][y][z]*A[0][x][y][z])
                - (lambda_ma*MA[0][x][y][z]*A[0][x][y][z])
                - m_A * A[0][x][y][z]
                + (d_a * laplacian(A,x,y,z))
                  ) * deltaT + A[0][x][y][z];

              if(A[1][x][y][z] != A[1][x][y][z]) 
              {
                //cout << "A\t(NaN)-> i:"<<i<<"\t-> ("<< x << y << z << ")" << A[1][x][y][z]<< "\n";
      //  return 1;
              }

            //Macrophages     
            
      /*****************************************************************
      * Assuming: contact with blood vessels only on one border bv = 0
      *           homogeneous contact with blood vessels bv = 1   
      *           contact with blood vessels given by function bv = 2
      ******************************************************************/
              source_mr = 0;            
              if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z))))    
                source_mr = alpha_mr * (m_estrela - MR[0][x][y][z]);
      /*****************************************************************/             
              MR[1][x][y][z] = ((- m_Mr * MR[0][x][y][z])
                     - (gamma_ma *MR[0][x][y][z] * A[0][x][y][z])
                     + (d_mr * laplacian(MR,x,y,z))
                     + source_mr ) * deltaT + MR[0][x][y][z];
                    
     
              if(MR[1][x][y][z] != MR[1][x][y][z] ) 
              {
                //cout << "MR\t(NaN)-> i: " << i << "-> "<< MR[1][x][y][z] <<"(" << x << y << z << ")" << "\n";
       // return 1;
              }
               MA[1][x][y][z] = ((-m_Ma * MA[0][x][y][z])
                     + (gamma_ma *MR[0][x][y][z] * A[0][x][y][z])
                     + (d_ma * laplacian(MA,x,y,z))
                      ) * deltaT + MA[0][x][y][z];
              if(MA[1][x][y][z] != MA[1][x][y][z] ) 
              {
                //cout << "MA\t(NaN)-> i: " << i << "-> "<< MA[1][x][y][z] <<"(" << x << y << z << ")"<< "\n" ;
      //  return 1;
              }

//*****************************************************************************

    //Simulates complete model
            }
            else if(simCase==0)
            // Add the equations here
            {
    //Antigenos
              A[1][x][y][z] = ( beta_A*A[0][x][y][z]*(1-(A[0][x][y][z]/k_A))
                    - ( lambda_mr*MR[0][x][y][z]*A[0][x][y][z])
                    - ( lambda_ma * MA[0][x][y][z] * A[0][x][y][z])
                    - ( lambda_afma*F[0][x][y][z]*A[0][x][y][z]*MA[0][x][y][z])
                    - ( lambda_afmr*F[0][x][y][z]*A[0][x][y][z]*MR[0][x][y][z])
                    - m_A * A[0][x][y][z]
                    + (d_a * laplacian(A,x,y,z))
                      ) * deltaT + A[0][x][y][z];
           
              if(A[1][x][y][z] < tol) 
              {
                A[1][x][y][z] = 0.0;
              }
              if(A[1][x][y][z] != A[1][x][y][z]) 
              {
                //cout << "A\t(NaN)-> i:"<<i<<" -> ("<< x << y << z << ") ->" << A[1][x][y][z] << "\n";
     //   return 1;
              }

          //Macrophages
            
      /*****************************************************************
      * Assuming: contact with blood vessels only on one border bv = 0
      *           homogeneous contact with blood vessels bv = 1   
      *           contact with blood vessels given by function bv = 2
      ******************************************************************/
              source_mr = 0;            
              if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z))))    
                source_mr = alpha_mr * (m_estrela - MR[0][x][y][z]);
      /*****************************************************************/  
               MR[1][x][y][z] = ((- m_Mr * MR[0][x][y][z])
                 - (gamma_ma * 
                 (k_m + k_MTNF * pow(TNF_alpha[0][x][y][z], h_MTNF)/(pow(n_MTNF, h_MTNF) +pow(TNF_alpha[0][x][y][z], h_MTNF)))
                        *pow(n_M10, h_M10)/(pow(n_M10, h_M10) +pow(IL_10[0][x][y][z], h_M10)) *
                  MR[0][x][y][z] * A[0][x][y][z])       
                 + (d_mr * laplacian(MR,x,y,z))
                 + source_mr ) * deltaT + MR[0][x][y][z];

              //if(MR[1][x][y][z] != MR[1][x][y][z] ) 
              //{
                //cout << "MR\t(NaN)-> i: " << i << "-> "<< MR[1][x][y][z] <<"(" << x << y << z << ")" << "\n";
     // return 1;
             // }

      
      /*****************************************************************
      * Assuming: contact with lymph vessels only on one border bv = 0
      *           homogeneous contact with lymph vessels bv = 1   
      *           contact with lymph vessels given by function bv = 2
      ******************************************************************/
              migration_ma = 0;            
              if (((lnv==0)&&(x==0))||(lnv==1)||((lnv==2)&&(is_lnvase(x,y,z))))    
                migration_ma = alpha_Ma * (MA_T - MA_L);
      /*****************************************************************/ 
              MA[1][x][y][z] = ((-m_Ma * MA[0][x][y][z])
                    + (gamma_ma * 
                      (k_m + k_MTNF * pow(TNF_alpha[0][x][y][z], h_MTNF)/(pow(n_MTNF, h_MTNF) +pow(TNF_alpha[0][x][y][z], h_MTNF)))
                        *pow(n_M10, h_M10)/(pow(n_M10, h_M10) +pow(IL_10[0][x][y][z], h_M10)) *
                      MR[0][x][y][z] * A[0][x][y][z])
                    + (d_ma * laplacian(MA,x,y,z))
                    - migration_ma ) * deltaT + MA[0][x][y][z];

          
              //if(MA[1][x][y][z] != MA[1][x][y][z] ) 
              //{
                //cout << "MA\t(NaN)-> i: " << i << "-> "<< MA[1][x][y][z] <<" -> (" << x << y << z << ")" << "\n";
    //   return 1;
              //}
     //Antibody
      
      /*****************************************************************
      * Assuming: contact with lymph vessels only on one border bv = 0
      *           homogeneous contact with lymph vessels bv = 1   
      *           contact with lymph vessels given by function bv = 2
      ******************************************************************/
              migration_f = 0;            
      //if (((lnv==0)&&(x==0))||(lnv==1)||((lnv==2)&&(is_lnvase(x,y,z))))    
              if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z))))
                migration_f = (alpha_f * (F_T - F_L));
      /*****************************************************************/ 
      
              F[1][x][y][z] = (
                   - ( lambda_afma * F[0][x][y][z] * A[0][x][y][z]*MA[0][x][y][z])
                   - ( lambda_afmr*F[0][x][y][z]*A[0][x][y][z]*MR[0][x][y][z])
                   - migration_f + (d_f * laplacian(F,x,y,z))
                     )* deltaT + F[0][x][y][z];

             // if(F[1][x][y][z] != F[1][x][y][z] ) 
             // {
              //  cout << "F\t(NaN)-> i: " << i << "-> "<< F[1][x][y][z] <<"(" << x << y << z << ")"<< "\n" ;
      //  return 1;
 
             // }
/********************************************************************/
              /*cytokines*/
              
               //TNF_alpha
               TNF_alpha[1][x][y][z] = (k_TNFM * pow(n_TNF6, h_TNF6)/(pow(n_TNF6, h_TNF6) +pow(IL_6[0][x][y][z], h_TNF6))
                * pow(n_TNF10,h_TNF10)/(pow(n_TNF10, h_TNF10) +pow(IL_10[0][x][y][z], h_TNF10)) * MA[0][x][y][z]
                - k_TNF * (TNF_alpha[0][x][y][z] - (MA[0][x][y][z]/MA[0][0][0][0]) * (MR[0][x][y][z] / MR[0][0][0][0]) * q_TNF) + d_TNF_alpha * laplacian(TNF_alpha,x,y,z)) * deltaT + TNF_alpha[0][x][y][z];
              
                //if(TNF_alpha[1][x][y][z] != TNF_alpha[1][x][y][z]) 
                //{
                  //cout << "TNF_alpha\t(NaN)-> i:"<<i<<" -> ("<< x << y << z << ") ->" << TNF_alpha[1][x][y][z] << "\n";
                   //   return 1;
               // }
                //IL_6
                 IL_6[1][x][y][z] = ((k_6m + k_6TNF * pow(TNF_alpha[0][x][y][z], h_6TNF)/(pow(n_6TNF, h_6TNF) +pow(TNF_alpha[0][x][y][z], h_6TNF)))
                                    * pow(n_66, h_66)/(pow(n_66, h_66) +pow(IL_6[0][x][y][z], h_66))
                                     * pow(n_610, h_610)/(pow(n_610, h_610) +pow(IL_10[0][x][y][z], n_610)) * MA[0][x][y][z]
                                    - k_6 * (IL_6[0][x][y][z] - (MA[0][x][y][z]/MA[0][0][0][0]) * (MR[0][x][y][z] / MR[0][0][0][0]) * q_IL6)) * deltaT + IL_6[0][x][y][z];
                // if(IL_6[1][x][y][z] != IL_6[1][x][y][z]) 
                //{
                 // cout << "IL_6\t(NaN)-> i:"<<i<<" -> ("<< x << y << z << ") ->" << IL_6[1][x][y][z] << "\n";
                    //   return 1;
                //}
                //IL_8
                 IL_8[1][x][y][z] = (k_8m + k_8TNF * pow(TNF_alpha[0][x][y][z], h_8TNF)/(pow(TNF_alpha[0][x][y][z], h_8TNF) +pow(n_8TNF, h_8TNF))
                  * pow(n_810, h_810)/(pow(n_810, h_810) +pow(IL_10[0][x][y][z], h_810)) * MA[0][x][y][z]
                                    - k_8 * (IL_8[0][x][y][z] - (MA[0][x][y][z]/MA[0][0][0][0]) * (MR[0][x][y][z] / MR[0][0][0][0]) * q_IL8)
                          ) * deltaT + IL_8[0][x][y][z]; 
               // if(IL_8[1][x][y][z] != IL_8[1][x][y][z]) 
               // {
               //   cout << "IL_8\t(NaN)-> i:"<<i<<" -> ("<< x << y << z << ") ->" << IL_8[1][x][y][z] << "\n";
                  //   return 1;
               // }
                //IL_10
                IL_10[1][x][y][z] = ((k_10m + k_106 * pow(IL_6[0][x][y][z], h_106)/(pow(n_106, h_106) +pow(IL_6[0][x][y][z], h_106))) * MA[0][x][y][z] 
                                       - (k_10 * (IL_10[0][x][y][z] - (MA[0][x][y][z]/MA[0][0][0][0]) * (MR[0][x][y][z] / MR[0][0][0][0]) * q_IL10)) 
                                     ) * deltaT + IL_10[0][x][y][z];
                //if(IL_10[1][x][y][z] != IL_10[1][x][y][z]) 
                //{
                  //cout << "IL_10\t(NaN)-> i:"<<i<<" -> ("<< x << y << z << ") ->" << IL_10[1][x][y][z] << "\n";
                   //   return 1;
                //}
            }
          }
        }
      }

      update(A);
      if(simCase==2||simCase==0)
      {
        update(MR);
        update(MA);
      }
      if(simCase==0)
      {
        update(F);
        update(TNF_alpha);
        update(IL_6);
        update(IL_8);
        update(IL_10);
      }
    }
    t++;
    //cout << "Limite: " << iterPerDay*days;
  }while(t < (iterPerDay*days));

  //if (saveFiles)
  //{
    //cout << "Closing files \n";
   // fclose(datamatlabT);
   // cout << "T closed \n";
   // fclose(datamatlabB);
   // cout << "B closed \n";
   // fclose(datamatlabP);
   // cout << "P closed \n";
   // fclose(datamatlabL);
   // cout << "L closed \n";
  //}

  cout << "\n" << Footer(t);
  double* VOI;
  VOI = new double [4];
  VOI[0] = TNF_alpha_T;
  VOI[1] = IL_6_T;
  VOI[2] = IL_8_T;
  VOI[3] = IL_10_T;
  //cout << "TNFalpha: " << TNF_alpha_T << "\n";
  //cout << "IL6: "  << IL_6_T << "\n";
  //cout << "IL8: "  << IL_8_T << "\n";
  //cout << "IL10: " <<  IL_10_T << "\n";
  return VOI;
}
