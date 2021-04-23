#include "QCircuit.h"
#include "tensorProd.h"
#include "matMul.h"


QCircuit::QCircuit(int nQReg) {
  nQubits = nQReg;
  nGates = 0;
  gateList.gates = {}; gateList.supports = {}; gateList.angles = {};
  isUComputed = false;
}

QCircuit::~QCircuit() {}

void QCircuit::addU1(int trg, double beta, double gamma, double delta) {
  gateList.gates.push_back("U1"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[3] = {beta, gamma, delta}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addPhase(int trg, double alpha)  {
  gateList.gates.push_back("Phase"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {alpha}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addH(int trg)  {
  gateList.gates.push_back("H"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addS(int trg) {
  gateList.gates.push_back("S"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addSdag(int trg) {
  gateList.gates.push_back("Sdag"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addT(int trg) {
  gateList.gates.push_back("T"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addTdag(int trg) {
  gateList.gates.push_back("Tdag"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addX(int trg) {
  gateList.gates.push_back("X"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addY(int trg){
  gateList.gates.push_back("Y"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addZ(int trg){
  gateList.gates.push_back("Z"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addRX(int trg, double theta) {
  gateList.gates.push_back("RX"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {theta}; 
  gateList.angles.push_back(angles); 
}

void QCircuit::addRY(int trg, double theta) {
  gateList.gates.push_back("RY"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {theta}; 
  gateList.angles.push_back(angles); 
}


void QCircuit::addRZ(int trg, double theta) {
  gateList.gates.push_back("RZ"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {theta}; 
  gateList.angles.push_back(angles); 
}


void QCircuit::addCNOT(int ctr, int trg) {
  gateList.gates.push_back("CNOT"); 
  int supports[2] = {ctr, trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
}


void QCircuit::compMatrix() {
}
void QCircuit::applyGates(Complex* state){
}
double QCircuit::expSigma(int trg, int sigma){
}
  

