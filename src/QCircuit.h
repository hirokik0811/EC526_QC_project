#ifndef QCIRCUIT_H
#define QCIRCUIT_H

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <bitset>  

#include "sparseMat.h"
#include "tensorProd.h"
#include "matMul.h"


using namespace std;


#define PI 3.141592653589793
#define  I  Complex(0.0, 1.0)
#define  Nfft 16
typedef complex <double> Complex;
//typedef struct{Complex* val; int* row; int* col; int nVal;} sparseCSR; 
typedef struct{vector<string> ops; vector<int*> supports; vector<double*> angles;} opList;


// For converting the gate name to integer allowing one to use switch for them
constexpr 
unsigned int str2int(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (str2int(str, h+1)*33) ^ str[h];
}

class QCircuit {
private:
  int nQubits, nGates;
  opList gateList; 
  sparseCSR U; // matrix representation of the circuit
  bool isUComputed; // to avoid computing the matrix twice
public: 
  QCircuit(int nQReg) {
    if (nQReg >= 32) { cout << "QCircuit support only <32 qubits. " << endl; return;}
    nQubits = nQReg;
    nGates = 0;
    isUComputed = false;
  }

  ~QCircuit() {}


  // Add a general U1 gate = RZ(beta)*RY(gamma)*RZ(delta)
  void addU1(int trg, double beta, double gamma, double delta) {
    if (trg >= nQubits) { cout << "The target qubit index " << trg << " is out of range. Please set it 0-"<< nQubits-1 << endl; return; }
    gateList.ops.push_back("U1"); 
    int* supports = new int[1]; supports[0] = trg; 
    gateList.supports.push_back(supports); 
    double* angles = new double[3]; angles[0] = beta; angles[1] = gamma; angles[2] = delta; 
    gateList.angles.push_back(angles); 
    nGates++; 
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
  }
  
  // Compute the matrix representation of one gate
  sparseCSR*  _gateMat(string* opName, int* supports, double* angles) {
    int d = (int)pow(2, nQubits); // dimension of the matrix rep
    sparseCSR gate; 
    // 2-qubit gate cases
    if (*opName == "CNOT") {
      int ctr = supports[0], trg = supports[1]; 
      Complex* val = new Complex[d]; // the values of the CNOT matrix
      int* row = new int[d+1]; 
      int* col = new int[d]; 
      for (int i = 0; i < d; ++i) {
        val[i] = 1.0; row[i] = i; 
        string bini = bitset<32>(i).to_string(); 
        if (bini[ctr] == '1' && bini[trg] == '0') { col[i] = i + (int)pow(2, trg); }
        else if (bini[ctr] == '1' && bini[trg] == '1') { col[i] = i-(int)pow(2, trg); }
        else { col[i] = i; }
      }
      row[d] = d; 
      gate.val = val;  gate.row = row;  gate.col = col; 
      gate.nVal = d; gate.nRows = d; gate.nCols = d; 
    }
    // 1-qubit gate cases
    else {
      int trg = supports[0];

      // identities supporting other qubits
      int firstIdDim = (int)pow(2, trg-1), secondIdDim = (int)pow(2, nQubits-trg-1); 
      sparseCSR firstId(firstIdDim, firstIdDim), secondId(secondIdDim, secondIdDim); 

      for (int i = 0; i < firstIdDim; ++i) { firstId.val[i] = 1.0; firstId.col[i] = i; firstId.row[i] = i; } 
      for (int i = 0; i < secondIdDim; ++i) { secondId.val[i] = 1.0; secondId.col[i] = i; secondId.row[i] = i; } 
      firstId.row[firstIdDim] = firstIdDim; secondId.row[secondIdDim] = secondIdDim; 

      // the target 1-qubit gate
      sparseCSR oneGate; 
      Complex* val; 
      int* row; 
      int* col;

      switch (str2int((*opName).c_str())) {
      case str2int("U1") :
        val = new Complex[4]; 
        row = new int[3]; 
        col = new int[4]; 
        val[0] = exp(I*(-angles[0]/2.0-angles[1]/2.0))*cos(angles[2]/2.0);
        val[1] = -exp(I*(-angles[0]/2.0+angles[1]/2.0))*sin(angles[2]/2.0);
        val[2] = exp(I*(angles[0]/2.0-angles[1]/2.0))*sin(angles[2]/2.0);
        val[3] = exp(I*(angles[0]/2.0+angles[1]/2.0))*cos(angles[2]/2.0);
        col[0] = 0; col[1] = 1; col[2] = 0; col[3] = 1; 
        row[0] = 0; row[1] = 2; row[2] = 4; 
        oneGate.nVal = 4;
      case str2int("Phase") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = exp(I*angles[0]);
        val[1] = exp(I*angles[0]);
        col[0] = 0; col[1] = 1;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      case str2int("H") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = 1/sqrt(2.0);
        val[1] = 1/sqrt(2.0);
        val[2] = 1/sqrt(2.0);
        val[3] = -1/sqrt(2.0);
        col[0] = 0; col[1] = 1; col[2] = 0; col[3] = 1;
        row[0] = 0; row[1] = 2; row[2] = 4; 
        oneGate.nVal = 4;
      case str2int("S") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = 1.0;
        val[1] = I;
        col[0] = 0; col[1] = 1;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      case str2int("Sdag") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = 1.0;
        val[1] = -I;
        col[0] = 0; col[1] = 1;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      case str2int("T") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = 1.0;
        val[1] = exp(I*PI/4.0);
        col[0] = 0; col[1] = 1;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      case str2int("Tdag") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = 1.0;
        val[1] = exp(-I*PI/4.0);
        col[0] = 0; col[1] = 1;
        row[0] = 0; row[1] = 1; row[2] = 2; 
      case str2int("X") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = 1.0;
        val[1] = 1.0;
        col[0] = 1; col[1] = 0;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      case str2int("Y") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = -I;
        val[1] = I;
        col[0] = 1; col[1] = 0;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      case str2int("Z") : 
        val = new Complex[2]; 
        row = new int[3]; 
        col = new int[2]; 
        val[0] = 1.0;
        val[1] = -1.0;
        col[0] = 0; col[1] = 1;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      case str2int("RX") :
        val = new Complex[4]; 
        row = new int[3]; 
        col = new int[4]; 
        val[0] = cos(angles[0]/2.0);
        val[1] = -I*sin(angles[0]/2.0);
        val[2] = cos(angles[0]/2.0);
        val[3] = -I*sin(angles[0]/2.0);
        col[0] = 0; col[1] = 1; col[2] = 0; col[3] = 1; 
        row[0] = 0; row[1] = 2; row[2] = 4; 
        oneGate.nVal = 2;
      case str2int("RY") :
        val = new Complex[4]; 
        row = new int[3]; 
        col = new int[4]; 
        val[0] = cos(angles[0]/2.0);
        val[1] = sin(angles[0]/2.0);
        val[2] = sin(angles[0]/2.0);
        val[3] = cos(angles[0]/2.0);
        col[0] = 0; col[1] = 1; col[2] = 0; col[3] = 1; 
        row[0] = 0; row[1] = 2; row[2] = 4; 
        oneGate.nVal = 2;
      case str2int("RZ") : 
        val = new Complex[4]; 
        row = new int[3]; 
        col = new int[4]; 
        val[0] = 1.0;
        val[1] = exp(I*angles[0]/2.0);
        col[0] = 0; col[1] = 1;
        row[0] = 0; row[1] = 1; row[2] = 2; 
        oneGate.nVal = 2;
      }
      oneGate.val = val;  oneGate.row = row;  oneGate.col = col; 
      oneGate.nRows = 2;  oneGate.nCols = 2; 
      
      if (trg != 0) { 
        gate = *tensorProdSparse(&oneGate, &firstId); 
      }  else {
        sparseCSR temp(oneGate); 
        gate = temp; 
      } 
      if (trg != nQubits-1) { 
        sparseCSR prevRes = gate; 
        gate = *tensorProdSparse(&secondId, &gate); 
        prevRes.freeMemory(); 
      }
      oneGate.freeMemory();  firstId.freeMemory();  secondId.freeMemory();
     
    }
    return &gate; 
  }

  // Compute the matrix representation of the circuit. 
  void compMatrix() {
    int d = (int)pow(2, nQubits); // dimension of the matrix rep
    sparseCSR matRep(); 
    for (int i = 0; i < nGates; ++i) {
      sparseCSR* gate = new sparseCSR;
      gate = _gateMat(&(gateList.ops[i]), gateList.supports[i], gateList.angles[i]); 

      // WRITE CODE FOR MULTIPLICATION 

    }
  }

  // Apply gates to the state vector
  void applyGates(Complex* state){
  }
  
  // Return the number of the qubits. 
  int getNQubits(){
    return nQubits;
  }
  
  int getNGates() { 
    return nGates; 
  }
  
  // Return the list of the gates. 
  opList getGateList(){
    return gateList;
  }

};



  
#endif 

