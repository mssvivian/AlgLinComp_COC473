#include <bits/stdc++.h>
#include <omp.h>
#include "trabalho1.cpp"

vector<vector<double>> generateConvergentMatrix(int n)
{
    // Cria a matriz nxn inicializada com 0
    vector<vector<double>> matrix(n, vector<double>(n, 0));
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        matrix[i][i] = n-1; // Define o valor na diagonal principal como 4
        if (i > 0)
        {
            matrix[i][i - 1] = -1; // Define o valor na diagonal inferior como -1
        }
        if (i < n - 1)
        {
            matrix[i][i + 1] = -1; // Define o valor na diagonal superior como -1
        }
    }

    return matrix;
}

int main()
{
    int size = 5;

    // vector<vector<double>> A = {
    //     {1.000, 2.000},
    //     {3.000, 4}};
    // vector<double> b = {0.1, 0.1};

    vector<vector<double>> A = generateConvergentMatrix(size);
    vector<double> b(size);
    b[0] = 1;
    

    convergencia(A);
    vector<double> x(b.size());
    double start, end;

    cout << "LU decomposition:" << endl;
    start = omp_get_wtime();
    x = luDecompositionSolve(A, b);
    end = omp_get_wtime();
    cout << "x: ";
    for (int i = 0; i < b.size(); i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << "Tempo de execução: " << end - start << " segundos." << endl;
    cout << endl;
    cout << endl;

    cout << "Gaussian Elimination:" << endl;
    start = omp_get_wtime();
    x = gaussianElimination(A, b);
    end = omp_get_wtime();
    cout << "x: ";
    for (int i = 0; i < b.size(); i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << "Tempo de execução: " << end - start << " segundos." << endl;
    cout << endl;
    cout << endl;

    cout << "Gauss-Seidel:" << endl;
    start = omp_get_wtime();
    x = gaussSeidel(A, b, 20, 1e-3, true);
    end = omp_get_wtime();
    cout << "x: ";
    for (int i = 0; i < b.size(); i++) {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << "Tempo de execução: " << end - start << " segundos." << endl;
    cout << endl;
    cout << endl;

    cout << "Jacobi:" << endl;
    start = omp_get_wtime();
    x = jacobiSolve(A, b, 20, 1e-3, true);
    end = omp_get_wtime();
    cout << "x: ";
    for (int i = 0; i < b.size(); i++) {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << "Tempo de execução: " << end - start << " segundos." << endl;


    // graph shit

    // cout << sizeof(double) << endl;
    // cout << "Jacobi:" << endl;
    // cout << fixed << setprecision(5);
    // cout << setw(20) << "tamanho" << setw(20) << " tempo" << endl;
    // for (int i = 1; i<=15; i++){
    //     vector<vector<double>> A = generateConvergentMatrix(1<<i);
    //     vector<double> b(1<<i);
    //     vector<double> x(1<<i);
    //     b[0] = 1;
    //     start = omp_get_wtime();
    //     x = jacobiSolve(A, b, 100, 1e-3, false);
    //     end = omp_get_wtime();
    //     cout << setw(20) << (1<<i) << " " << setw(20) << end - start << endl;
    // }

    // cout << "Gauss-Sidel:" << endl;
    // cout << fixed << setprecision(5);
    // cout << setw(20) << "tamanho" << setw(20) << " tempo" << endl;
    // for (int i = 1; i<=15; i++){
    //     vector<vector<double>> A = generateConvergentMatrix(1<<i);
    //     vector<double> b(1<<i);
    //     vector<double> x(1<<i);
    //     b[0] = 1;
    //     start = omp_get_wtime();
    //     x = gaussSeidel(A, b, 100, 1e-3, false);
    //     end = omp_get_wtime();
    //     cout << setw(20) << (1<<i) << " " << setw(20) << end - start << endl;
    // }

    // cout << "Gauss:" << endl;
    // cout << fixed << setprecision(5);
    // cout << setw(20) << "tamanho" << setw(20) << " tempo" << endl;
    // for (int i = 1; i<=12; i++){
    //     vector<vector<double>> A = generateConvergentMatrix(1<<i);
    //     vector<double> b(1<<i);
    //     vector<double> x(1<<i);
    //     b[0] = 1;
    //     start = omp_get_wtime();
    //     x = gaussianElimination(A, b);
    //     end = omp_get_wtime();
    //     cout << setw(20) << (1<<i) << " " << setw(20) << end - start << endl;
    // }
    
    // cout << "LU:" << endl;
    // cout << fixed << setprecision(5);
    // cout << setw(20) << "tamanho" << setw(20) << " tempo" << endl;
    // for (int i = 1; i<=12; i++){
    //     vector<vector<double>> A = generateConvergentMatrix(1<<i);
    //     vector<double> b(1<<i);
    //     vector<double> x(1<<i);
    //     b[0] = 1;
    //     start = omp_get_wtime();
    //     x = luDecompositionSolve(A, b);
    //     end = omp_get_wtime();
    //     cout << setw(20) << (1<<i) << " " << setw(20) << end - start << endl;
    // }
    return 0;
}
