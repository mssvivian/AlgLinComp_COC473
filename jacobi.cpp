#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

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

int main(){
    vector<vector<double>> A = {
        {3, -1, -1},
        {-1, 3, -1},
        {-1, -1, 3}
    };
    vector<double> b = {1.0, 2.0, 1.0};
    vector<double> x(b.size());
    auto start = omp_get_wtime();
    x = jacobiSolve(A, b);
    auto end = omp_get_wtime();
    cout << "Tempo de execução: " << end - start << " segundos." << endl;
    return 0;
}