#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

using namespace std;

void DecomposicaoLU(vector<vector<double>> A, vector<double> b){
    int n = A.size();

    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<vector<double>> U(n, vector<double>(n, 0.0));
    vector<double> y(n, 0.0);
    vector<double> x(n, 0.0);

    // Decompondo A em L e U
    for (int i = 0; i < n; i++){
        for (int j = i; j < n; j++){
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++){
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
        for (int j = i +1; j < n; j++){
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++){
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
        L[i][i] = 1;
    }

    // Ly = b
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }

    // Ux = y
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    cout << "Solução x:\n";
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i << "] = " << x[i] << "\n";
    }

}

int main() {
    vector<vector<double>> A = {
        {1, 2, 2},
        {4, 4, 2},
        {4, 6, 4}
    };

    vector<double> b = {3, 6, 10};

    DecomposicaoLU(A, b);

    return 0;
}
