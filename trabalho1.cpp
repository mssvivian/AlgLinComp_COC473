#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

bool convergencia(const vector<vector<double>> &A)
{
    int N = A.size();
    for (int i = 0; i < N; i++)
    {
        double count_row = 0.0, count_column = 0.0;
        for (int j = 0; j < N; j++)
        {
            if (j != i)
            {
                count_row += fabs(A[i][j]);
                count_column += fabs(A[j][i]);
            }
        }
        if (fabs(A[i][i]) < count_row && fabs(A[i][i]) < count_column)
        {
            cout << "Nao diagonalmente dominante.\n";
            return false;
        }
    }
    cout << "Diagonalmente dominante.\n";
    return true;
}

vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b)
{
    int N = A.size();
    vector<double> x(N);

    for (int i = 0; i < N; i++)
    {
        if (A[i][i] == 0)
        {
            int swap_row = -1;
            for (int j = i + 1; j < N; j++)
            {
                if (A[j][i] != 0)
                {
                    swap_row = j;
                    break;
                }
            }
            if (swap_row == -1)
            {
                cout << "Sistema sem solucao\n";
                return {};
            }
            swap(A[i], A[swap_row]);
            swap(b[i], b[swap_row]);
        }

        double pivot = A[i][i];

        #pragma omp parallel for
        for (int j = i + 1; j < N; j++)
        {
            double factor = A[j][i] / pivot;
            #pragma omp parallel for
            for (int k = i; k < N; k++)
            {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    for (int i = N - 1; i >= 0; i--)
    {
        double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
        for (int j = i + 1; j < N; j++)
        {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

vector<double> luDecompositionSolve(vector<vector<double>> A, vector<double> b)
{
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    vector<double> y(n, 0), x(n, 0);

    for (int k = 0; k < n - 1; k++)
    {
#pragma omp parallel for
        for (int i = k + 1; i < n; i++)
        {
            A[i][k] = A[i][k] / A[k][k];
        }
#pragma omp parallel for
        for (int j = k + 1; j < n; j++)
        {
            for (int i = k + 1; i < n; i++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }

    // Ly = b
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
        for (int j = 0; j < i; j++)
        {
            sum += A[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }

    // Ux = y
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
        for (int j = i + 1; j < n; j++)
        {
            sum += A[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / A[i][i];
    }

    return x;
}

vector<double> gaussSeidel(vector<vector<double>> &A, vector<double> &b, int max_iter = 20, double tol = 1e-3, bool print = false)
{
    int n = A.size();
    vector<double> x(n, 1.0); // Inicializa x com 1.0
    vector<double> x_last(n, 1.0);
    vector<double> x_diff(n);

    // cout << fixed << setprecision(5);
    if(print) cout << setw(20) << "Iteração" << setw(20) << "Residuo" << endl;

    for (int iter = 0; iter < max_iter; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    sum += A[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        // Check convergence
        double diff_norm = 0.0;
#pragma omp parallel for reduction(+ : diff_norm)
        for (int i = 0; i < n; i++)
        {
            diff_norm += (x_last[i] - x[i]) * (x_last[i] - x[i]);
        }
        diff_norm = sqrt(diff_norm);

        double norm = 0.0;
#pragma omp parallel for reduction(+ : norm)
        for (int i = 0; i < n; i++)
        {
            norm += (x[i]) * (x[i]);
        }
        norm = sqrt(norm);

        double residue = (diff_norm / norm);

        if (print)
        {
            cout << setw(10) << iter + 1;

            // Print residue
            cout << setw(20) << residue << endl;
        }

        if (residue < tol)
        {
            if(print){
                cout << "Convergiu em " << iter + 1 << " iterações." << endl;
                cout << "Residuo: " << residue << endl;
            }
            break;
        }

        x_last = x;
    }

    return x;
}

vector<double> jacobiSolve(vector<vector<double>> &A, vector<double> &b, int max_iter = 20, double tol = 1e-3, bool print = false)
{
    int n = A.size();
    vector<double> x(n, 1.0); // Inicializa x com 1.0
    vector<double> x_new(n);
    vector<double> x_diff(n);

    // cout << fixed << setprecision(5);
    if(print) cout << setw(20) << "Iteração" << setw(20) << "Residuo" << endl;

    for (int iter = 0; iter < max_iter; iter++)
    {
#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        // Check convergence
        double diff_norm = 0.0;
#pragma omp parallel for reduction(+ : diff_norm)
        for (int i = 0; i < n; i++)
        {
            diff_norm += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }
        diff_norm = sqrt(diff_norm);

        double norm = 0.0;
#pragma omp parallel for reduction(+ : norm)
        for (int i = 0; i < n; i++)
        {
            norm += (x_new[i]) * (x_new[i]);
        }
        norm = sqrt(norm);

        double residue = (diff_norm / norm);

        if (print)
        {
            cout << setw(10) << iter + 1;

            // Print residue
            cout << setw(20) << residue << endl;
        }

        if (residue < tol)
        {
            if(print){
                cout << "Convergiu em " << iter + 1 << " iterações." << endl;
                cout << "Residuo: " << residue << endl;
            }
            break;
        }

        x = x_new;
    }

    return x;
}
