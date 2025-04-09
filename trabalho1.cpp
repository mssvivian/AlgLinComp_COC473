#include <iostream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>

using namespace std;

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

vector<double> luDecompositionSolve(vector<vector<double>> A, vector<double> b) {
    int N = A.size();
    vector<vector<double>> L(N, vector<double>(N, 0));
    vector<vector<double>> U(N, vector<double>(N, 0));
    vector<double> y(N, 0), x(N, 0);

    
    // Decompondo A em L e U
    for (int i = 0; i < N; i++){
        for (int j = i; j < N; j++){
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++){
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
        for (int j = i +1; j < N; j++){
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++){
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
        L[i][i] = 1;
    }

    // Ly = b
    for (int i = 0; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }

    // Ux = y
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

vector<double> gaussSeidel(vector<vector<double>> A, vector<double> B, vector<double> X, double TOL, int iteracoes_max) {
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
            cout << k << " iteracoes\n";
            return X;
        }

        Xo = X;
        k++;
    }

    cout << "Numero maximo de iteracoes atingido.\n";
    return X;
}

vector<double> jacobiSolve(vector<vector<double>> &A, vector<double> &b, int max_iter = 20, double tol = 1e-3) {
    int n = A.size();
    vector<double> x(n, 1.0); // Inicializa x com 1.0
    vector<double> x_new(n);
    vector<double> x_diff(n);

    cout << fixed << setprecision(5);
    cout << setw(20) << "Iteração" << setw(20) << "x anterior" << setw(20) << "x novo" << setw(20) << "Residuo" << endl;

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

        cout << setw(10) << iter + 1;

        // Print x (previous)
        cout << setw(20);
        if (n <= 6) {
            for (int i = 0; i < n; i++) {
                cout << x[i] << " ";
            }
        } else {
            for (int i = 0; i < 2; i++) {
                cout << x[i] << " ";
            }
            cout << "... ";
            for (int i = n - 1; i < n; i++) {
                cout << x[i] << " ";
            }
        }

        // Print x_new (current)
        cout << setw(20);
        if (n <= 6) {
            for (int i = 0; i < n; i++) {
                cout << x_new[i] << " ";
            }
        } else {
            for (int i = 0; i < 2; i++) {
                cout << x_new[i] << " ";
            }
            cout << "... ";
            for (int i = n - 1; i < n; i++) {
                cout << x_new[i] << " ";
            }
        }

        // Print residue
        cout << setw(20) << residue << endl;

        if (residue < tol) {
            cout << "Convergiu em " << iter + 1 << " iterações." << endl;
            cout << "Residuo: " << residue << endl;
            cout << "Solução: ";
            for (int i = 0; i < n; i++) {
                cout << x_new[i] << " ";
            }
            break;
        }

        x.swap(x_new);
    }

    return x;
}

int main() {
    double TOL;
    int iteracoes_max;
    int N;
    cout << "Digite o tamanho da matriz: ";
    cin >> N;

    vector<vector<double>> A(N, vector<double>(N));
    vector<double> B(N), X(N, 1);

    cout << "Digite a matriz A (" << N << "x" << N << "):\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Digite o vetor B:\n";
    for (int i = 0; i < N; i++) {
        cin >> B[i];
    }


    int escolha;
    cout << "Escolha o metodo:\n";
    cout << "1 - Gauss-Seidel (Iterativo)\n";
    cout << "2 - Eliminacao de Gauss (Direto)\n";
    cout << "3 - Fatoracao LU (Direto)\n";
    cout << "4 - Jacobi (Iterativo)\n";
    cin >> escolha;

    vector<double> resultado;
    if (escolha == 1) {
        cout << "Digite a tolerancia: ";
        cin >> TOL;
        cout << "Digite o numero maximo de iteracoes: ";
        cin >> iteracoes_max;
        if (convergencia(A)) {
            resultado = gaussSeidel(A, B, X, TOL, iteracoes_max);
        } else {
            cout << "A matriz nao e diagonalmente dominante, o metodo pode nao convergir.\n";
        }
    } else if (escolha == 2) {
        resultado = gaussianElimination(A, B);
    } else if (escolha == 3) {
        resultado = luDecompositionSolve(A, B);
    } else if (escolha == 4){
        cout << "Digite a tolerancia: ";
        cin >> TOL;
        cout << "Digite o numero maximo de iteracoes: ";
        cin >> iteracoes_max;
        if (convergencia(A)) {
            resultado = jacobiSolve(A, B, TOL, iteracoes_max);
        } else {
            cout << "A matriz nao e diagonalmente dominante, o metodo pode nao convergir.\n";
        }
    } else {
        cout << "Opcao invalida.\n";
        return 1;
    }

    if (!resultado.empty()) {
        cout << "Vetor X (solucao):\n";
        for (double xi : resultado) {
            cout << xi << " ";
        }
        cout << '\n';
    }

    return 0;
}

//vector<vector<double>> A = { {3, -1, -1}, {-1, 3, -1}, {-1, -1, 3} };
//vector<double> B = {1, 2, 1};
//vector<double> X = {1, 1, 1};
// 1,25 1,5 1,25
