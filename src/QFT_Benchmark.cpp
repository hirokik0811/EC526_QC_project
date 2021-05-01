/* Benchmarking N-qubit quantum fourier transforms */

#include "QCircuit.h"
#include <iomanip>
#include <fstream>
#include <chrono>

void cRz(QCircuit* circ, int ctr, int trg, double thet) {
    circ->addRZ(trg, thet / 2.0);
    circ->addCNOT(ctr, trg);
    circ->addRZ(ctr, thet / 2.0);
    circ->addRZ(trg, -thet / 2.0);
    circ->addCNOT(ctr, trg);
}

void swap(QCircuit* circ, int ctr, int trg) {
    circ->addCNOT(trg, ctr);
    circ->addCNOT(ctr, trg);
    circ->addCNOT(trg, ctr);
}


int main() {

    ofstream outfile;
    outfile.open("Benchmark.txt");

    bool printQCircuitData = false;
    bool compExact = false;
    bool printOutputs = false;

    std::chrono::time_point<std::chrono::steady_clock> begin_time;
    std::chrono::time_point<std::chrono::steady_clock> end_time;
    std::chrono::duration<double> difference_in_time;
    double difference_in_seconds;

    outfile << std::fixed;
    outfile << setprecision(9);

    int N2;
    for (int N = 2; N <= 10; N++)
    {
        N2 = pow(2, N);
        /* CONSTRUCT FT QCIRCUIT */
        QCircuit circ(N);

        int phasePower = 1;
        for (int i = N - 1; i >= 0; i--)
        {
            circ.addH(i);
            phasePower = 1;
            for (int j = i - 1; j >= 0; j--)
            {
                cRz(&circ, j, i, PI / pow(2.0, phasePower));
                phasePower++;
            }
        }

        for (int i = 0; i < N / 2; i++)
        {
            swap(&circ, i, N - 1 - i);
        }


        /* PRINT QCIRCUIT DATA */
        if (printQCircuitData)
        {
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
        }


        /* COMPUTE THE EXACT MATRIX REPRESENTATION */
        if (compExact)
        {
            MatrixXC ftExact(N2, N2);
            Complex omega = exp(I * 2.0 * PI / (double)N2);
            for (int i = 0; i < N2; ++i) {
                for (int j = 0; j < N2; ++j) {
                    ftExact(i, j) = pow(omega, i * j) / sqrt(N2);
                }
            }
            cout << endl;
            cout << "difference with the exact matrix" << endl;
            cout << circ.getUDense() - ftExact << endl;
        }

        /* COMPUTE THE MATRIX REPRESENTATION OF FT CIRCUIT */
        cout << "\nComputing matrix representation of a " << N << "-Qubit Quantum Fourier Transform..." << endl;
        outfile << N;

        begin_time = std::chrono::steady_clock::now();
        circ.getUSparse();
        end_time = std::chrono::steady_clock::now();
        difference_in_time = end_time - begin_time;
        difference_in_seconds = difference_in_time.count();

        cout << "Completed in " << difference_in_seconds << " seconds." << endl;
        outfile << "\t" << difference_in_seconds;

        SpMat state(N2, 1);
        vector<T> tripletList;
        tripletList.push_back(T(0, 0, 1.0));
        state.setFromTriplets(tripletList.begin(), tripletList.end());

        begin_time = std::chrono::steady_clock::now();
        circ.outputState(state);
        end_time = std::chrono::steady_clock::now();
        difference_in_time = end_time - begin_time;
        difference_in_seconds = difference_in_time.count();

        cout << "Output computed in " <<  difference_in_seconds << " seconds given an input of |0>." << endl;
        outfile << "\t" << difference_in_seconds << endl;


        /* COMPUTE THE OUTPUT STATE FOR EACH INPUT STATE */
        if (printOutputs)
        {
            for (int i = 0; i < N2; i++) {
                tripletList[0] = T(i, 0, 1.0);
                state.setFromTriplets(tripletList.begin(), tripletList.end());
                cout << "Output vector with |" << i << "> input: " << endl;
                cout << MatrixXC(circ.outputState(state)) << endl;
            }
        }
    }

    return 0;
}


