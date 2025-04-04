#include <stdio.h>
#include <math.h>

#define N 3        
#define TOL 0.001   
#define iteracoes_max  1000 

int convergencia(double A[N][N]) { 
    int i, j;
    double count_row, count_column;
    
    for (i = 0; i < N; i++) {
        count_row = 0.0;
        count_column = 0.0;
        for (j = 0; j < N; j++) {
            if (j != i) {
                count_row += fabs(A[i][j]);
                count_column += fabs(A[j][i]);
            }
        }
        if (fabs(A[i][i]) < count_row && fabs(A[i][i]) < count_column) {
            printf("Nao convergente.\n");
            return 0; 
        }
    }
    printf("Diagonalmente dominante.\n");
    return 1;
}

void gaussSeidel(double A[N][N], double B[N], double X[N]) {
    double Xo[N];  
    int i, j;
    int k = 1;
    double sigma;

    for (i = 0; i < N; i++) {
        Xo[i] = X[i];
    }

    while (k <= iteracoes_max ) {
        for (i = 0; i < N; i++) {
            sigma = 0.0;
            for (j = 0; j < i; j++) {
                sigma += A[i][j] * X[j];
            }
            for (j = i + 1; j < N; j++) {
                sigma += A[i][j] * Xo[j]; 
            }
            X[i] = (B[i] - sigma) / A[i][i];
        }
        int conver = 1;
        for (i = 0; i < N; i++) {
            if (fabs(X[i] - Xo[i]) > TOL) {
                conver = 0;
                break;
            }
        }
        if (conver) {
            printf("%d iteracoes\n", k);
            return;
        }
        for (i = 0; i < N; i++) {
            Xo[i] = X[i];
        }
        k = k + 1; 
    }

    printf("Numero maximo de iteracoes atingido.\n");
}

int main() {
    double A[N][N] = { {3, -1, -1}, {-1, 3, -1}, {-1, -1, 3} };
    double B[N] = {1, 2, 1};
    double X[N] = {1, 1, 1}; 

    if (convergencia(A)) { 
        gaussSeidel(A, B, X);
    } else {
        printf("A matriz não é diagonalmente dominante. Métodos iterativos podem não convergir.\n");
    }

    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.4f\n", i + 1, X[i]);
    }

    return 0;
}