#include "IS_Model.h"
#include "fstream"
#define N 51

using namespace std;

int main(){
  //IS_Model m(1,1);
//  IS_Model* m = new IS_Model();
  IS_Model* model = new IS_Model();
  //model->setSaveFiles(1);
  //model->setSimulationCase(1);

  double *parameters = new double[N];
  //valores iniciais sao todos iguais a 0?

  double a0 = 2.0;
  parameters[0] = a0;

  //Verificar quais destes devem ser removidos
  double TNF_0 =  0.04;
  parameters[1] = TNF_0;
  double IL_10_0 = 0.00027;
  parameters[2] = IL_10_0;
  double IL_6_0 = 0.04;
  parameters[3] = IL_6_0;
  double IL_8_0 = 0.02;
  parameters[4] = IL_8_0;

  double n_610 = 34.8/1000; //pg/mm^3
  parameters[5] = n_610;
  double n_66 = 560/1000; // pg/mm^3
  parameters[6] = n_66;
  double n_6TNF = 185/1000; // pg/mm^3
  parameters[7] = n_6TNF;
  double n_810 = 17.4/1000; // pg/mm63
  parameters[8] = n_810;
  double n_8TNF = 185/1000; // pg/mm^3
  parameters[9] = n_8TNF;
  double n_106 = 560/1000; // pg/mm^3
  parameters[10] = n_106;
  double n_TNF10 = 17.4/1000; // pg/mm^3
  parameters[11] = n_TNF10;
  double n_TNF6 = 560/1000; // pg/mm^3
  parameters[12] = n_TNF6;
  double n_M10 = 4.35/1000; // pg/mm^3
  parameters[13] = n_M10;
  double n_MTNF = 0.1; // pg/mm^3
  parameters[14] = n_MTNF;

  double h_106 = 3.68;
  parameters[15] = h_106;
  double h_6TNF = 2;
  parameters[16] = h_6TNF;
  double h_66 = 1;
  parameters[17] = h_66;
  double h_610 = 4;
  parameters[18] = h_610;
  double h_8TNF = 3;
  parameters[19] = h_8TNF;
  double h_810 = 1.5;
  parameters[20] = h_810;
  double h_TNF10 = 3;
  parameters[21] = h_TNF10;
  double h_TNF6 = 2;
  parameters[22] = h_TNF6;
  double h_M10 = 0.3;
  parameters[23] = h_M10;
  double h_MTNF = 3.16;
  parameters[24] = h_MTNF;

  double A_T     = a0/1000;
  parameters[25] = A_T;
  
  double TNF_alpha_T = TNF_0;
  parameters[26] = TNF_alpha_T;
  double IL_6_T = IL_6_0;
  parameters[27] = IL_6_T;
  double IL_8_T = IL_8_0;
  parameters[28] = IL_6_T;
  double IL_10_T = IL_10_0;
  parameters[29] = IL_10_T;
  double d_a        = 0.00037;       // antigen diffusion (Haessler) mm^3/day
  parameters[30] = d_a;
  double d_mr       = 0.0432;        //resting macrophage diffusion (estimated) mm^3/day
  parameters[31] = d_mr;
  double d_ma       = 0.3;           //macrophage diffusion (estimated) mm^3/day
  parameters[32] = d_ma;
  double d_f        = 0.016;          //antibody diffusion (estimated) mm^3/day
  parameters[33] = d_f;
//TNF
  double k_TNFM = 0.6 * 24/(1000); //pg/mm^3*day* # of cells
  parameters[34] = k_TNFM;
  double k_TNF = 24; // day^-1
  parameters[35] = k_TNF;
  double q_TNF = 4; //pg/mm^3
  parameters[36] = q_TNF;
//IL_6
  double k_6m = 0.81 * 24/1000; //pg/mm^3*day* # of cells
  parameters[37] = k_6m;
  double k_6TNF = 0.81 * 24 /1000; //pg/mm^3*day* # of cells
  parameters[38] = k_6TNF;
  double k_6 = 4.64*24;//day^-1
  parameters[39] = k_6;
  double q_IL6 = 4; //pg/mm^3
  parameters[40] = q_IL6;
//IL_8
  double k_8m = 0.56 *24/1000;//pg/mm^3*day* # of cells
  parameters[41] = k_8m;
  double k_8TNF = 0.56*24/1000; //pg/mm^3*day* # of cells
  parameters[42] = k_8TNF;
  double k_8 = 4.64*24; //day^-1
  parameters[43] = k_8;
  double q_IL8 = 2; //pg/mm^3
  parameters[44] = q_IL8;
//IL_10
  double k_10m = 0.0191*24/1000; //pg/mm^3*day* # of cells
  parameters[45] = k_10m;
  double k_106 = 0.0191 *24/1000; //pg/mm^3*day* # of cells
  parameters[46] = k_106;
  double k_10 = 0.5*24; //day^-1
  parameters[47] = k_10;
  double q_IL10 = 0.27; // pg/mm^3
  parameters[48] = q_IL10;

  double k_m = 0.0414 *24; // 1/day
  parameters[49] = k_m;
  double k_MTNF = 8.65 * 24; // 1/day
  parameters[50] = k_MTNF;

  double* param_ = new double [4];//Grava temporariamente o valor para o paramentro*90
  double* paramM = new double [4];//Grava temporariamente o valor para o paramentro*110
  double paramOriginal;//Grava temporariamente o valor original do parametro
  double* indiceSensibilidade = new double [4];//Guarda os índices de sensibilidade de cada variável
  double* resultOriginal = new double [4];
  // chamar solve com parametros base
  resultOriginal = model->solve(parameters);
  // Arquivo de saida
  fstream arquivo;
  arquivo.open("output.txt");
  //TNF_alpha_T, IL_6_T, IL_8_T, IL_10_T ordem do retorno do solve
  for(int i = 0; i< N; i++){
    paramOriginal = parameters[i];
    parameters[i] = paramOriginal*0.90;
    param_ = model->solve(parameters);//solve -10% na variável i
    parameters[i] = paramOriginal*1.10;
    paramM = model->solve(parameters);//solve +10% na variável i
    parameters[i] = paramOriginal;
    for(int k = 0; k<4; k++){
      //Conta para calcular o indice de sensibilidade
      //Barbara:inseri esse if para evitar divisoes por zero
      if((paramM[k] - param_[k])==0){
        indiceSensibilidade[k] = 0.0;
      }else{// Barbara: eu tinha digitado uma divisao no lugar de multiplicacao
        indiceSensibilidade[k] = ((paramM[k] - param_[k])/(2*0.1*paramOriginal))*
                              (paramOriginal/resultOriginal[k]);    
      }
    }
    //resultOriginal ESTÁ RETORNANDO 0
    arquivo << i << ", " << indiceSensibilidade[0] << ", " << indiceSensibilidade[1] 
    << ", " << indiceSensibilidade[2] << ", " << indiceSensibilidade[3] << endl;
    paramOriginal = 0.0;
  }
  //free
  delete[] resultOriginal; 
  return 0;
}
