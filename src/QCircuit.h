#ifndef QCIRCUIT_H
#define QCIRCUIT_H

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <bitset>  

//#include "sparseMat.h"
#include "tensorProd.h"
#include "matMul.h"


using namespace std;


#define PI 3.141592653589793
#define  I  Complex(0.0, 1.0)
#define  Nfft 16
typedef complex <double> Complex;
typedef struct{vector<string> ops; vector<int*> supports; vector<double*> angles;} opList;

typedef Matrix<Complex, Dynamic, Dynamic> MatrixXC; // the complex dense matrix for printing and debugging 


// For converting the gate name to integer allowing one to use switch for them
constexpr 
unsigned int str2int(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (str2int(str, h+1)*33) ^ str[h];
}

class QCircuit {
private:
  /************** member variables ***************/ 
  int nQubits, nGates;
  opList gateList; 
  int idxComputed; 
  int d; // Full dimension of the matrix representation
  SpMat U; // matrix representation of the circuit
  bool isUComputed; // to avoid computing the matrix twice
  
  /************** member functions ***************/ 
  // Construct an identity matrix with dim x dim dimensions
  SpMat _identityMat(int dim) {
    vector<T> tripletList; // the non-zero values
    for (int i = 0; i < dim; ++i) { tripletList.push_back(T(i, i, 1.0)); }
    SpMat id(dim, dim); 
    id.setFromTriplets(tripletList.begin(), tripletList.end());
    return id; 
  }
  
  // Compute the matrix representation of one gate
  SpMat _gateMat(string* opName, int* supports, double* angles) {
    SpMat gate(d, d); 
    // 2-qubit gate cases
    if (*opName == "CNOT") {
      int ctr = supports[0], trg = supports[1]; 
      vector<T> tripletList; // the non-zero values of the CNOT matrix
      for (int i = 0; i < d; ++i) {
        string bini = bitset<32>(i).to_string(); 
        if (bini[31-ctr] == '1' && bini[31-trg] == '0') { tripletList.push_back(T(i, i + (int)pow(2, trg), 1.0)); }
        else if (bini[31-ctr] == '1' && bini[31-trg] == '1') { tripletList.push_back(T(i, i-(int)pow(2, trg), 1.0)); }
        else { tripletList.push_back(T(i, i, 1.0)); }
      }
      gate.setFromTriplets(tripletList.begin(), tripletList.end());
    }
    // 1-qubit gate cases
    else {
      int trg = supports[0];

      // identities supporting other qubits
      int firstIdDim = (int)pow(2, trg), secondIdDim = (int)pow(2, nQubits-trg-1); 
      SpMat firstId = _identityMat(firstIdDim), secondId = _identityMat(secondIdDim); 
      

      // the target 1-qubit gate
      vector<T> tripletList; // the non-zero values of the 1-qubit gate

      switch (str2int((*opName).c_str())) {
      case str2int("U1") :
        tripletList.push_back(T(0, 0, exp(I*(-angles[0]/2.0-angles[1]/2.0))*cos(angles[2]/2.0)));
        tripletList.push_back(T(0, 1, -exp(I*(-angles[0]/2.0+angles[1]/2.0))*sin(angles[2]/2.0)));
        tripletList.push_back(T(1, 0, exp(I*(angles[0]/2.0-angles[1]/2.0))*sin(angles[2]/2.0)));
        tripletList.push_back(T(1, 1, exp(I*(angles[0]/2.0+angles[1]/2.0))*cos(angles[2]/2.0)));
        break;
      case str2int("Phase") :  
        tripletList.push_back(T(0, 0, exp(I*angles[0])));
        tripletList.push_back(T(1, 1, exp(I*angles[0])));
        break;
      case str2int("H") : 
        tripletList.push_back(T(0, 0, 1.0/sqrt(2.0)));
        tripletList.push_back(T(0, 1, 1.0/sqrt(2.0)));
        tripletList.push_back(T(1, 0, 1.0/sqrt(2.0)));
        tripletList.push_back(T(1, 1, -1.0/sqrt(2.0)));
        break;
      case str2int("S") : 
        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, I));
        break;
      case str2int("Sdag") : 
        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, -I));
        break;
      case str2int("T") : 
        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, exp(I*PI/4.0)));
        break;
      case str2int("Tdag") : 
        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, exp(-I*PI/4.0)));
        break;
      case str2int("X") : 
        tripletList.push_back(T(0, 1, 1.0));
        tripletList.push_back(T(1, 0, 1.0));
        break;
      case str2int("Y") : 
        tripletList.push_back(T(0, 1, -I));
        tripletList.push_back(T(1, 0, I));
        break;
      case str2int("Z") : 
        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, -1.0));
        break;
      case str2int("RX") :
        tripletList.push_back(T(0, 0, cos(angles[0]/2.0)));
        tripletList.push_back(T(0, 1, -I*sin(angles[0]/2.0)));
        tripletList.push_back(T(1, 0, -I*sin(angles[0]/2.0)));
        tripletList.push_back(T(1, 1, cos(angles[0]/2.0)));
        break;
      case str2int("RY") :
        tripletList.push_back(T(0, 0, cos(angles[0]/2.0)));
        tripletList.push_back(T(0, 1, -sin(angles[0]/2.0)));
        tripletList.push_back(T(1, 0, sin(angles[0]/2.0)));
        tripletList.push_back(T(1, 1, cos(angles[0]/2.0)));
        break;
      case str2int("RZ") : 
        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, exp(I*angles[0])));
        break;
      }
      
      SpMat oneGate(2, 2);
      oneGate.setFromTriplets(tripletList.begin(), tripletList.end());

      SpMat right; 
      if (trg != 0) { 
        right = tensorProdSparse(oneGate, firstId);
      }  else {
        right = oneGate; 
      }
      if (trg != nQubits-1) { 
        gate = tensorProdSparse(secondId, right); 
      } else {
        gate = right; 
      }
    }
    return gate; 
  }

  // Compute the matrix representation of the circuit. 
  void _compMatrix() {
    if (!isUComputed) {
      int d = (int)pow(2, nQubits); // dimension of the matrix rep 
      for (int i = idxComputed+1; i < nGates; ++i) {
        SpMat gate;
        gate = _gateMat(&(gateList.ops[i]), gateList.supports[i], gateList.angles[i]); 
        
        if (i == 0) { U = gate; } 
        else { U = matrixMul(gate, U); }

      }
      idxComputed = nGates-1; 
      isUComputed = true; 
    }
  }

public: 
  QCircuit(int nQReg) {
    if (nQReg >= 32) { cout << "QCircuit support only <32 qubits. " << endl; return;}
    nQubits = nQReg;
    nGates = 0;
    idxComputed = -1; 
    d = (int)pow(2, nQubits); 
    U = _identityMat(d); 
    isUComputed = true;
  }

  ~QCircuit() {
    for (int i = 0; i < nGates; ++i) {
      delete[] gateList.supports[i]; delete[] gateList.angles[i]; 
    }
  }

  // Add a general U1 gate = RZ(beta)*RY(gamma)*RZ(delta)
  void addU1(int trg, double beta, double gamma, double delta) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("U1"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[3]; angles[0] = beta; angles[1] = gamma; angles[2] = delta; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add a trivial phase gate = exp(i*alpha)
  void addPhase(int trg, double alpha)  {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("Phase"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[1]; angles[0] = alpha; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add a Hadamard gate = 1/sqrt{2}*{{1, 1}, {1, -1}}
  void addH(int trg)  {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("H"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }


  // Add an S gate = {{1, 0}, {0, i}} = RZ(pi/2)
  void addS(int trg) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("S"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add an Sdag gate = {{1, 0}, {0, -i}} = RZ(-pi/2)
  void addSdag(int trg) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("Sdag"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }


  // Add a T gate = {{1, 0}, {0, exp(i*pi/4)}} = RZ(pi/4)
  void addT(int trg) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("T"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add a Tdag gate = {{1, 0}, {0, exp(-i*pi/4)}} = RZ(-pi/4)
  void addTdag(int trg) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("Tdag"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add an X gate = {{0, 1}, {1, 0}} = sigma_1
  void addX(int trg) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("X"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add a Y gate = {{0, -i}, {i, 0}} = sigma_2
  void addY(int trg){
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("Y"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add a Z gate = {{1, 0}, {0, -1}} = sigma_3
  void addZ(int trg){
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("Z"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add an RX gate = {{cos(theta/2), -i*sin(theta/2)}, {-i*sin(theta/2), cos(theta/2)}}
  void addRX(int trg, double theta) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("RX"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[1]; angles[0] = theta;
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add an RY gate = {{cos(theta/2), -sin(theta/2)}, {sin(theta/2), cos(theta/2)}}
  void addRY(int trg, double theta) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("RY"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[1]; angles[0] = theta;
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Add an RZ gate = {{1, 0}, {0, exp(i*theta)}}
  void addRZ(int trg, double theta) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("RZ"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[1]; angles[0] = theta;
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }


  // Add a CNOT gate
  void addCNOT(int ctr, int trg) {
    if (ctr >= nQubits) { cout << "The control qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("CNOT"); 
    int* supports = new int[1]; supports[0] = ctr; supports[1] = trg; 
    gateList.supports.push_back(supports);
    double* angles = new double[0]; 
    gateList.angles.push_back(angles); 
    nGates++; 
    isUComputed = false;
  }

  // Apply gates to the state vector
  SpMat outputState(SpMat inputState){
    _compMatrix();
    return matrixMul(U, inputState); 
  }
  
  // Compute the 1-qubit expectation value <in|U^dag Obs_i U|in>
  // allowed observables are sigmas: obs = "X", "Y", "Z"
  double expVal(string* obs, int trg, SpMat inputState) {
    if (obs[0] != "X" && obs[0] != "Y" && obs[0] != "Z") { 
       cout << "The allowed observables are only X, Y, or Z" << endl; 
       return -1.0; 
    }
    int * sup = new int[1]; sup[0] = trg; 
    double * ang = new double[0];
    SpMat observable = _gateMat(obs, sup, ang); 
    delete[] sup; delete[] ang; 
    SpMat outSt = outputState(inputState); 
    SpMat outStDag = outSt.conjugate().transpose(); 
    SpMat firstProd = matrixMul(outStDag, observable); 
    SpMat innProd = matrixMul(firstProd, outSt); 
    return (innProd.coeff(0, 0)).real(); 
  }

  // Simulate the measurement results of the computational basis with finite #samplings. Store the sampling result to countsPtr. 
  void ZBasisSampling(int* countsPtr, int nSamples, SpMat inputState) {
    SpMat outSt = outputState(inputState); 
    double* dist = new double[d]; 
    for (int i = 0; i < d; ++i) { 
      dist[i] = (outSt.coeff(i, 0)*conj(outSt.coeff(i, 0))).real();
    }
    
    for (int i = 0; i < d; ++i) { 
      countsPtr[i] = 0; 
    }

    // For random # generation
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);

    for (int i = 0; i < nSamples; ++i) { 
      double smp = distribution(generator); 
      for (int j = 0; j < d; ++j) { 
         smp -= dist[j]; 
         if (smp <= 0.0) {
           countsPtr[j]++; 
           break; 
         }
      }
    }
    delete[] dist; 
  }
  
  
  // Return the number of the qubits. 
  int getNQubits(){
    return nQubits;
  }
  
  // Return the number of the gates. 
  int getNGates() { 
    return nGates; 
  }
  
  // Return the list of the gates. 
  opList getGateList(){
    return gateList;
  }

  // Return the CSR sparse matrix representation of U
  SpMat getUSparse() { 
    _compMatrix(); 
    return U;
  } 

  // Return the dense matrix representation of U
  MatrixXC getUDense() { 
    _compMatrix();
    return MatrixXC(U); 
  }

};



  
#endif 

