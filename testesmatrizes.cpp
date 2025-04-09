#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <time.h>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>


using namespace std;

vector<vector<double>> generateMatrix(double n) {
    // Cria a matriz nxn inicializada com 0
    vector<vector<double>> matrix(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i) {
        matrix[i][i] = 4; // Define o valor na diagonal principal como 4
        if (i > 0) {
            matrix[i][i - 1] = -1; // Define o valor na diagonal inferior como -1
        }
        if (i < n - 1) {
            matrix[i][i + 1] = -1; // Define o valor na diagonal superior como -1
        }
    }

    return matrix;
}

vector<double> createVector(int size, int position) {
    // Cria o vector<double> inicializado com zeros
    vector<double> vec(size, 0.0);

    // Verifica se a posição está dentro dos limites do vetor
    if (position >= 0 && position < size) {
        vec[position] = -1; // Define -1 na posição especificada
    } else {
        cerr << "Erro: posição fora dos limites do vetor!" << endl;
    }

    return vec;
}

bool convergencia(const vector<vector<double>>& A) { 
    int N = A.size();
    for (int i = 0; i < N; i++) {
        double count_row = 0.0, count_column = 0.0;
        for (int j = 0; j < N; j++) {
            if (j != i) {
                count_row += fabs(A[i][j]);
                count_column += fabs(A[j][i]);
            }
        }
        if (fabs(A[i][i]) < count_row && fabs(A[i][i]) < count_column) {
            cout << "Nao convergente.\n";
            return false; 
        }
    }
    cout << "Diagonalmente dominante.\n";
    return true;
}

vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int N = A.size();
    vector<double> x(N);

    for (int i = 0; i < N; i++) {
        if (A[i][i] == 0) {
            int swap_row = -1;
            for (int j = i + 1; j < N; j++) {
                if (A[j][i] != 0) {
                    swap_row = j;
                    break;
                }
            }
            if (swap_row == -1) {
                cout << "Sistema sem solucao\n";
                return {};
            }
            swap(A[i], A[swap_row]);
            swap(b[i], b[swap_row]);
        }

        double pivot = A[i][i];
        for (int j = i; j < N; j++) {
            A[i][j] /= pivot;
        }
        b[i] /= pivot;

        for (int j = i + 1; j < N; j++) {
            double fator = A[j][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= fator * A[i][k];
            }
            b[j] -= fator * b[i];
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }

    return x;
}

int gaussSeidel(vector<vector<double>> A, vector<double> B, vector<double> X, double TOL, int iteracoes_max) {
    int N = A.size();
    vector<double> Xo = X;
    int k = 1;

    while (k <= iteracoes_max) {
        for (int i = 0; i < N; i++) {
            double sigma = 0.0;
            for (int j = 0; j < i; j++) {
                sigma += A[i][j] * X[j];
            }
            for (int j = i + 1; j < N; j++) {
                sigma += A[i][j] * Xo[j];
            }
            X[i] = (B[i] - sigma) / A[i][i];
        }

        bool conver = true;
        for (int i = 0; i < N; i++) {
            if (fabs(X[i] - Xo[i]) > TOL) {
                conver = false;
                break;
            }
        }
        if (conver) {
            // cout << k << " iteracoes\n";
            return k ;
        }

        Xo = X;
        k++;
    }

    cout << "Numero maximo de iteracoes atingido.\n";
    return k;
}

void DecomposicaoLU(vector<vector<double>> A, vector<double> b){
    int n = A.size(); //show de bola esse size 

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

    // imprimindo vetor x
    // cout << "Solução x:\n";
    // for (int i = 0; i < n; i++) {
    //     std::cout << "x[" << i << "] = " << x[i] << "\n";
    // }

}

int jacobiSolve(vector<vector<double>> &A, vector<double> &b, int max_iter = 20, double tol = 1e-3) {
    int n = A.size();
    vector<double> x(n, 1.0); // Inicializa x com 1.0
    vector<double> x_new(n);
    vector<double> x_diff(n);

    // cout << fixed << setprecision(5);
    // cout << setw(20) << "Iteração" << setw(20) << "x anterior" << setw(20) << "x novo" << setw(20) << "Residuo" << endl;

    for (int iter = 0; iter < max_iter; iter++) {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            #pragma omp parallel for reduction(+:sum)
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        // Check convergence
        double diff_norm = 0.0;
        #pragma omp parallel for reduction(+:diff_norm)
        for (int i = 0; i < n; i++) {
            diff_norm += (x_new[i] - x[i])*(x_new[i] - x[i]);
        }
        diff_norm = sqrt(diff_norm);

        double norm = 0.0;
        #pragma omp parallel for reduction(+:norm)
        for (int i = 0; i < n; i++) {
            norm += (x_new[i])*(x_new[i]);
        }
        norm = sqrt(norm);

        double residue = (diff_norm / norm);

        // cout << setw(10) << iter + 1;

        // Print x (previous)
        // cout << setw(20);
        // if (n <= 6) {
        //     for (int i = 0; i < n; i++) {
        //         cout << x[i] << " ";
        //     }
        // } else {
        //     for (int i = 0; i < 2; i++) {
        //         cout << x[i] << " ";
        //     }
        //     cout << "... ";
        //     for (int i = n - 1; i < n; i++) {
        //         cout << x[i] << " ";
        //     }
        // }

        // Print x_new (current)
        // cout << setw(20);
        // if (n <= 6) {
        //     for (int i = 0; i < n; i++) {
        //         cout << x_new[i] << " ";
        //     }
        // } else {
        //     for (int i = 0; i < 2; i++) {
        //         cout << x_new[i] << " ";
        //     }
        //     cout << "... ";
        //     for (int i = n - 1; i < n; i++) {
        //         cout << x_new[i] << " ";
        //     }
        // }

        // Print residue
        // cout << setw(20) << residue << endl;

        if (residue < tol) {
            // cout << "Convergiu em " << iter + 1 << " iterações." << endl;
            // cout << "Residuo: " << residue << endl;
            // cout << "Solução: ";
            // for (int i = 0; i < n; i++) {
            //     cout << x_new[i] << " ";
            // }
            return iter + 1;
        }

        x.swap(x_new);
    }

    return 0;
}

int main() {

    // vector<vector<double>> B = generateMatrix(10000);

    // vector<double> b = createVector(10000,4999);

    // vector<double> X(10000, 1);

    vector<int> tamanho = {1024, 2048, 4096, 8192, 10000}; 

    for (int i = 0; i < 5; i++){
        int x = tamanho[i];

        vector<vector<double>> A = generateMatrix(x);
        vector<double> b = createVector(x,(x/2)-1);
        vector<double> X(x, 1);

        clock_t start_time, end_time;
        double cpu_time_used;
        start_time = clock();

        gaussianElimination(A, b);

        end_time = clock();
        cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
        printf("Tempo de execução: %10f segundos\n", cpu_time_used);
        // cout << k << "\n";     
    }



    return 0;
}


// {
//  {4, -1, 0,-1, 0},
//  {-1, 4, -1, 0, -1},
//  {0, -1, 4, -1, 0},
//  {0, 0, -1, 4, -1},
//  {0, 0, 0, -1, 4}
// };