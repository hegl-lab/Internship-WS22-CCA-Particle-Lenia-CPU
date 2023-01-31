#ifndef PARTICLE_LENIA_CPU_PARAMETERS_H
#define PARTICLE_LENIA_CPU_PARAMETERS_H


struct Parameters {
    double mu_k = 4.0;
    double sigma_k = 1.0;
    double w_k = 0.022;
    double mu_g = 0.6;
    double sigma_g = 0.15;
    double c_rep = 1.0;
    double dt = 0.1;
    double h = 0.01;
};


#endif //PARTICLE_LENIA_CPU_PARAMETERS_H
