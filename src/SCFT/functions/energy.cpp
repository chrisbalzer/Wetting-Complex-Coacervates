#include "../SCFT.h"

void SCFT::calcEnergy(double* phi_k, double& f) {
      double f_ID, f_DH, f_DHChain, f_EV;

      IdealEnergy(phi_k,f_ID);
      EVEnergy(phi_k,f_EV);
      DHEnergy(phi_k,f_DH);
      DHChainEnergy(phi_k,f_DHChain);

      f = f_ID + f_EV + f_DH + f_DHChain;
}
