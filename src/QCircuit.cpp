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


// Add a general U1 gate = RZ(beta)*RY(gamma)*RZ(delta)
void QCircuit::addU1(int trg, double beta, double gamma, double delta) {
  gateList.gates.push_back("U1"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[3] = {beta, gamma, delta}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add a trivial phase gate = exp(i*alpha)
void QCircuit::addPhase(int trg, double alpha)  {
  gateList.gates.push_back("Phase"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {alpha}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add a Hadamard gate = 1/sqrt{2}*{{1, 1}, {1, -1}}
void QCircuit::addH(int trg)  {
  gateList.gates.push_back("H"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}


// Add an S gate = {{1, 0}, {0, i}} = RZ(pi/2)
void QCircuit::addS(int trg) {
  gateList.gates.push_back("S"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add an Sdag gate = {{1, 0}, {0, -i}} = RZ(-pi/2)
void QCircuit::addSdag(int trg) {
  gateList.gates.push_back("Sdag"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}


// Add a T gate = {{1, 0}, {0, exp(i*pi/4)}} = RZ(pi/4)
void QCircuit::addT(int trg) {
  gateList.gates.push_back("T"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add a Tdag gate = {{1, 0}, {0, exp(-i*pi/4)}} = RZ(-pi/4)
void QCircuit::addTdag(int trg) {
  gateList.gates.push_back("Tdag"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add an X gate = {{0, 1}, {1, 0}} = sigma_1
void QCircuit::addX(int trg) {
  gateList.gates.push_back("X"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports); 
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add a Y gate = {{0, -i}, {i, 0}} = sigma_2
void QCircuit::addY(int trg){
  gateList.gates.push_back("Y"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add a Z gate = {{1, 0}, {0, -1}} = sigma_3
void QCircuit::addZ(int trg){
  gateList.gates.push_back("Z"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add an RX gate = {{cos(theta/2), -i*sin(theta/2)}, {-i*sin(theta/2), cos(theta/2)}}
void QCircuit::addRX(int trg, double theta) {
  gateList.gates.push_back("RX"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {theta}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add an RY gate = {{cos(theta/2), -sin(theta/2)}, {sin(theta/2), cos(theta/2)}}
void QCircuit::addRY(int trg, double theta) {
  gateList.gates.push_back("RY"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {theta}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Add an RZ gate = {{1, 0}, {0, exp(i*theta)}}
void QCircuit::addRZ(int trg, double theta) {
  gateList.gates.push_back("RZ"); 
  int supports[1] = {trg}; 
  gateList.supports.push_back(supports);
  double angles[1] = {theta}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}


// Add a CNOT gate
void QCircuit::addCNOT(int ctr, int trg) {
  gateList.gates.push_back("CNOT"); 
  int supports[2] = {ctr, trg}; 
  gateList.supports.push_back(supports);
  double angles[0] = {}; 
  gateList.angles.push_back(angles); 
  nGates++; 
}

// Compute the matrix representation of the circuit. 
void QCircuit::compMatrix() {
}

// Apply gates to the state vector
void QCircuit::applyGates(Complex* state){
}

  

