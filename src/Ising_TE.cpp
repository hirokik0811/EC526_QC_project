/* Demonstration of 4-qubit time evolution of <0|exp(iHt) Z_0 exp(-iHt) |0> where H is the four-site transverse Ising model with PBC. exp(-iHt) is approximated by the Trotter-Suzuki expansion.   */  

#include <fstream>
#include "QCircuit.h"

#define N_QUBITS 4


void RZZ(QCircuit* circ, int ctr, int trg, double thet) { 
  circ->addCNOT(ctr, trg);
  circ->addPhase(trg, -thet/2.0); 
  circ->addRZ(trg, thet);
  circ->addCNOT(ctr, trg);
}

void TrotterStep(QCircuit* circ, double deltaT, double J, double g) {
  for (int i = 0; i < N_QUBITS; ++i) {
     RZZ(circ, i, (i+1)%N_QUBITS, -2.0*J*deltaT); 
  }
  for (int i = 0; i < N_QUBITS; ++i) { 
     circ->addRX(i, -2.0*g*deltaT); 
  }
}


int main() { 
  FILE * File = fopen ("time_evolution_TFIM.txt","w");
  
  QCircuit circ(N_QUBITS);
  vector<T> tripletList; 
  tripletList.push_back(T(0, 0, 1.0));
  SpMat state((int)pow(2, N_QUBITS), 1); 
  state.setFromTriplets(tripletList.begin(), tripletList.end());
  string obs = "Z"; 
  
  double stopT = 10; 
  double measure_step = 0.5; 
  int iter_per_step = 5; 
  
  circ.expVal(&obs, 0, state); 
  fprintf(File, "{{0, %g} ", circ.expVal(&obs, 0, state));
  for (int i = 0; i < (int)(stopT/measure_step); ++i) {
    for (int j = 0; j < iter_per_step; ++j) { 
       TrotterStep(&circ, measure_step/(double)iter_per_step, 1.0, 1.0); 
    }
    fprintf(File, ", {%g, %g}", (double)(i+1)*measure_step, circ.expVal(&obs, 0, state)); 
  }
  fprintf(File, "} "); 
  fclose(File); 
  
  
  //cout << circ.getUDense() << endl; 
  return 0;
}
    

