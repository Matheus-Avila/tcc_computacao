#ifndef _IS_Model_H_
#define _IS_Model_H_

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

//each 10 space equals 1 cm (1 × 10^−2 m)
const int    Xspace   = 10.0; // 10 mm = 1 cm
const int    Yspace   = 10.0; // 10 mm = 1 cm
const int    Zspace   = 10.0; // 10 mm = 1 cm
const int    SPACE    = Xspace * Yspace * Zspace; // 1 cm^3
const int    IC_SPACE = 0.4*(Xspace * Yspace * Zspace); //initial condition
const int    source   = 100*pow(10,0);
const double SCALE    = pow(10,-3);
const double MOL      = 6.02*pow(10,23);
const int    buffer   = 2;

class IS_Model{

  private:

    double A [buffer][Xspace][Yspace][Zspace];
    double MR [buffer][Xspace][Yspace][Zspace];
    double MA [buffer][Xspace][Yspace][Zspace];
    double F [buffer][Xspace][Yspace][Zspace];

    //CREA basic function
    double IL_6 [buffer][Xspace][Yspace][Zspace];
    double TNF_alpha [buffer][Xspace][Yspace][Zspace];
    double IL_8 [buffer][Xspace][Yspace][Zspace];
    double IL_10 [buffer][Xspace][Yspace][Zspace];

    //CREA regulation function

    //MA&MR
    //double HU_M_TNF [buffer][Xspace][Yspace][Zspace];
    //double HD_M_IL10 [buffer][Xspace][Yspace][Zspace];

    //TNF
    //double HD_TNF_IL6 [buffer][Xspace][Yspace][Zspace];
    //double HD_TNF_IL10 [buffer][Xspace][Yspace][Zspace];

    //IL_6
    //double HU_IL6_TNF [buffer][Xspace][Yspace][Zspace];
    //double HD_IL6_IL6 [buffer][Xspace][Yspace][Zspace];
    //double HD_IL6_IL10 [buffer][Xspace][Yspace][Zspace];

    //IL_8
    //double HU_IL8_TNF [buffer][Xspace][Yspace][Zspace];
    //double HD_IL8_IL10 [buffer][Xspace][Yspace][Zspace];

    //IL_10
    //double HU_IL10_IL6 [buffer][Xspace][Yspace][Zspace];




    //CREA regulation function constants.
    double n_610;
    double n_66;
    double n_6TNF;
    double n_810;
    double n_8TNF;
    double n_106;
    double n_TNF10;
    double n_TNF6;
    double n_M10;
    double n_MTNF;

    double h_106;
    double h_6TNF;
    double h_66;
    double h_610;
    double h_8TNF;
    double h_810;
    double h_TNF10;
    double h_TNF6;
    double h_M10;
    double h_MTNF;

    int simCase;
    int days;
    int points;
    double deltaT;
    double iterPerDay;
    double deltaX, deltaY, deltaZ;

    int lnv;
    int bv;
    double tol;
    double source_mr;
    double migration_ma;
    double migration_f;
    /*
    double migration_TNF;
    double migration_IL_6;
    double migration_IL_8;
    double migration_IL_10;
    */


    double m0;
    double a0;
    double th0;
    double b0;
    double p0;
    double f0;

    //CREA
    double TNF_0;
    double IL_10_0;
    double IL_8_0;
    double IL_6_0;

    double t_estrela;
    double b_estrela;
    double p_estrela;
    double f_estrela;
    double m_estrela;
    // time
    double MA_T, MA_L, MR_T;
    double Th, B, P;
    double F_T, F_L, A_T;

    double IL_6_T, TNF_alpha_T, IL_8_T, IL_10_T;

    // difussion constant
    double d_a;
    double d_mr;
    double d_ma;

    double d_IL_6;
    double d_TNF_alpha;
    double d_IL_8;
    double d_IL_10;

    //the constant CREA
    // MR & MA
    double k_m;
    double k_MTNF;

    //IL6
    double k_6m;
    double k_6TNF;
    double k_6;
    double q_IL6;
    //TNF
    double k_TNFM;
    double k_TNF;
    double q_TNF;
    //IL8
    double k_8m;
    double k_8TNF;
    double k_8;
    double q_IL8;
    //IL10
    double k_10m;
    double k_106;
    double k_10;
    double q_IL10;
    // CREA constant ENDS

    double beta_A;
    double k_A;
    double m_A;
    double m_Mr;
    double m_Ma;
    double alpha_Ma;
    double gamma_ma;
    double lambda_mr;
    double lambda_ma;
    double lambda_afmr;
    double lambda_afma;
    double b_th;
    double b_p;
    double b_pb;
    double b_pp;
    double ro_t;
    double ro_b;
    double ro_p;
    double ro_f;
    double alpha_t;
    double alpha_b;
    double alpha_p;
    double alpha_f;
    double alpha_mr;
    double d_f;

    int saveFiles;
    char *dir;
    FILE* datamatlabA;
    FILE* datamatlabMr;
    FILE* datamatlabMa;
    FILE* datamatlabT;
    FILE* datamatlabB;
    FILE* datamatlabP;
    FILE* datamatlabF;
    FILE* datamatlabL;
    FILE* datamatlabTNF_alpha;
    FILE* datamatlabIL_8;
    FILE* datamatlabIL_10;
    FILE* datamatlabIL_6;
    FILE* datamatlabCytokine;
    std::string Header();
    std::string Footer(long int t);
    int checkFile(FILE* theFile);
    int calcIntegral(double vec[][Xspace][Yspace][Zspace], double *V);
    int calcIntegral_lv(double vec[][Xspace][Yspace][Zspace], double *V);
    int calcIntegral_bv(double vec[][Xspace][Yspace][Zspace], double *V);
    void initialize();
    void update(double vec[][Xspace][Yspace][Zspace]);
    double laplacian(double vec[][Xspace][Yspace][Zspace], int x, int y, int z);
    int is_bvase(int x, int y, int z);
    int is_lnvase(int x, int y, int z);

  public:
    IS_Model();
    //IS_Model(int sCase, int sFile);
    //~IS_Model();
    void setSaveFiles(int sf);
    void setSimulationCase(int sc);
    double* solve(double* parameters);

};

#endif
