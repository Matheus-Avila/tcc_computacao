#ifndef _solve_h_
#define _solve_h_

class solve{
    private:
    double *macrofagos;
    double *oligodendrocitos;
    double *citocinas;
    double deltaX, deltaT;
    int tempo_num_pt;

    public:
    solve();
    ~solve();
    double* resolve(double* parametros);
};

#endif