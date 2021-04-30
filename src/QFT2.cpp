/* Demonstration of 2-qubit quantum fourier transform */  


#include "QCircuit.h"

void cRz(QCircuit* circ, int ctr, int trg, double thet) {
  circ->addRZ(trg, thet/2.0); 
  circ->addCNOT(ctr, trg); 
  circ->addRZ(ctr, thet/2.0); 
  circ->addRZ(trg, -thet/2.0);
  circ->addCNOT(ctr, trg);
}


int main() {
  QCircuit circ(2);
  circ.addH(1); 
  cRz(&circ, 0, 1, PI/2.0); 
  circ.addH(0);
  circ.addCNOT(1, 0); 
  circ.addCNOT(0, 1); 
  circ.addCNOT(1, 0); 
  
  
  cout << "Gate List: " << endl;
  opList gateList = circ.getGateList(); 
  cout << "Length of the gate list : " << gateList.ops.size() << endl;
  cout << "Length of the support list : " << gateList.supports.size() << endl;
  cout << "Length of the angle list : " << gateList.angles.size() << endl;
  for (int it = 0; it < circ.getNGates(); ++it) {
    cout << "gate: " << gateList.ops[it] << endl;
    cout << "zeroth support: " << (gateList.supports[it])[0] << endl;
    cout << "zeroth angle: " << (gateList.angles[it])[0] << endl;
  }
  cout << "circuit matrix: " << endl;
  cout << circ.getUDense() << endl;
  
  MatrixXC ft4(4, 4); 
  Complex omega = exp(I*PI/2.0); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) {
      ft4(i, j) = pow(omega, i*j)/2.0; 
    }
  }
  cout << endl; 
  cout << "difference with the exact matrix" << endl; 
  cout << circ.getUDense() - ft4 << endl;

  
  vector<T> tripletList; 
  tripletList.push_back(T(0, 0, 1.0));
  SpMat state(4, 1); 
  state.setFromTriplets(tripletList.begin(), tripletList.end());
  cout << "Output vector with |0> input: " << endl;
  cout << MatrixXC(circ.outputState(state)) << endl;
  
  
  tripletList[0] = T(1, 0, 1.0);
  state.setFromTriplets(tripletList.begin(), tripletList.end());
  cout << "Output vector with |1> input: " << endl;
  cout << MatrixXC(circ.outputState(state)) << endl;
  
  tripletList[0] = T(2, 0, 1.0);
  state.setFromTriplets(tripletList.begin(), tripletList.end());
  cout << "Output vector with |2> input: " << endl;
  cout << MatrixXC(circ.outputState(state)) << endl;
  
  tripletList[0] = T(3, 0, 1.0);
  state.setFromTriplets(tripletList.begin(), tripletList.end());
  cout << "Output vector with |3> input: " << endl;
  cout << MatrixXC(circ.outputState(state)) << endl;
  return 0;
}
    

