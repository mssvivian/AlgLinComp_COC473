#include <iostream>
#include <Eigen/Dense>

// Fazer a troca de linhas caso o pivô seja zero (evitar divisão por zero)
void swapLines(Eigen::MatrixXd A, Eigen::VectorXd b,int n,int max){

  while(n<= max){
    if (A(n+1,n+1) != 0){
      Eigen::VectorXd temp = A.row(n+1);
      A.row(n+1) = A.row(n);
      A.row(n) = temp;
      double temp1 = b(n+1);
      b(n+1) = b(n);
      b(n) = temp1;
      break;
    }
    n++;
  }
}

int main() {
  // Define o tamanho do sistema linear (número de equações e variáveis)
  const int max = 3;

  // Cria a matriz de coeficientes A (n x n)
  Eigen::MatrixXd A(max, max);

  // Cria o vetor de termos independentes b (n x 1)
  Eigen::VectorXd b(max);

  // Inicializa a matriz A e o vetor b com valores
  A << 2, 1, 1,
       1, 3, 2,
       1, 0, 4;

  b << 4, 5, 6;
  int n = 0;

  while(n<max){
    double pivo = A(n,n);
    if(pivo == 0){
      swapLines(A,b,n,max);
      pivo = A(n,n); // pivô mudará em caso de swap
    }
    A.row(n) = A.row(n)/pivo;
    b(n) = b(n)/pivo;
    for (int i = n+1; i < max; i++) { 
      double fator = A(i, n);
      Eigen::VectorXd temp = A.row(n);
      A.row(i) -= temp*fator;
      b(i) -= b(n) * fator;
    }
    n++;
  }
  std::cout << "Matrix:\n" << A << std::endl;
  std::cout << "Vetor:\n" << b << std::endl;
  return 0;
}