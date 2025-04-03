#include <iostream>
#include <vector>
#include <cmath>

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