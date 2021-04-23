#ifndef QCIRCUIT_H
#define QCIRCUIT_H

#include <iostream>
#include <cmath>
#include <complex>
#include <list>
#include <string>
using namespace std;


#define PI 3.141592653589793
#define  I  Complex(0.0, 1.0)
#define  Nfft 16
typedef complex <double> Complex;
typedef struct{Complex* val, int* row, int* col, int nVal} sparseCSR; 
typedef struct{list<string> gates, list<int*> supports, list<double*> angles} opList;

class QCircuit {
  int nQubits, nGates;
  opList gateList; 
  sparseCSR U; // matrix representation of the circuit
  bool isUComputed; // to avoid computing the matrix twice
public: 
  QCircuit(int nQReg);
  ~QCircuit();
  void addU1(int trg, double beta, double gamma, double delta);
  void addPhase(int trg, double alpha); 
  void addH(int trg);
  void addS(int trg);
  void addSdag(int trg);
  void addT(int trg);
  void addTdag(int trg);
  void addX(int trg);
  void addY(int trg);
  void addZ(int trg);
  void addRX(int trg, double theta); 
  void addRY(int trg, double theta);
  void addRZ(int trg, double theta);
  void addCNOT(int ctr, int trg);
  void compMatrix();
  void applyGates(Complex* state);
}
  
#endif 

