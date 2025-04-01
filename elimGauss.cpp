
//#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
using namespace std;


int main() {
  // Define o tamanho do sistema linear (número de equações e variáveis)
  const int size = 3;

  // Cria a matriz de coeficientes A (n x n)
  // Eigen::MatrixXd A(size, size);
  vector<vector<double>> A(size, vector<double>(size));
  vector<double> b(size);
  vector<double> x(size);


  // Cria o vetor de termos independentes b (n x 1)
  // Eigen::VectorXd b(size);

  // Preenche a matriz A e o vetor b com valores aleatórios
  // srand(time(0));
  // #pragma omp parallel for
  // for (int i = 0; i < size; i++) {
  //   #pragma omp parallel for
  //   for (int j = 0; j < size; j++) {
  //     A[i][j] = rand() % 100 + 1;
  //   }
  //   b[i] = rand() % 100 + 1;
  // }

  A = { {1, 0, 1},
        {0, 1, 0},
        {1, 0, 2}};

  b = {1, 2, 3};

  //auto start = omp_get_wtime();

  for(int i = 0; i<size; i++){
    // Check if the pivot element is zero, if so, swap with a non-zero row
    if(A[i][i] == 0) {
      int swap_row = -1;
      for(int j = i + 1; j < size; j++) {
        if(A[j][i] != 0) {
          swap_row = j;
          break;
        }
      }
      if(swap_row == -1) {
        cout << "Sistema sem solução" << endl;
        return 0;
      }

      swap(A[i], A[swap_row]);
      swap(b[i], b[swap_row]);
    }

    double pivot = A[i][i];
    for(int j = i; j < size; j++) {
      A[i][j] /= pivot;
    }

    b[i] /= pivot;

    for(int j = i + 1; j < size; j++) {
      double fator = A[j][i];
      for(int k = i; k < size; k++) {
          A[j][k] -= fator * A[i][k];
      }
      b[j] -= fator * b[i];
    }

  }

  // Perform back substitution to find the solution vector x
  for (int i = size - 1; i >= 0; i--) {
    x[i] = b[i];
    for (int j = i + 1; j < size; j++) {
      x[i] -= A[i][j] * x[j];
    }
  }

  //auto end = omp_get_wtime();

  cout << "Matriz A:\n";
  for(int i = 0; i<size; i++){
    for(int j = 0; j<size; j++){
      cout << A[i][j] << " ";
    }
    cout << '\n';
  }
  cout << '\n';


  cout << "Vetor b:\n";
  for (int i = 0; i < size; i++) {
    cout << b[i] << " ";
  }
  cout << '\n';

  cout << "Vetor x (solução):\n";
  for (int i = 0; i < size; i++) {
    cout << x[i] << " ";
  }
  cout << '\n';

  //cout << "Tempo de execução: " << end-start << "s" << endl;

  // cout << "Matrix:\n" << A << endl;
  // cout << "Vetor:\n" << b << endl;
  return 0;
}