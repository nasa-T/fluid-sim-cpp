#include "fluid.h"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
// #include <eigen3/Eigen/src/Core/util/ForwardDeclarations.h>
#include <eigen3/Eigen/Eigen>

bool adjacent(int i1, int j1, int i2, int j2) {
    return (abs(i1 - i2) + abs(j1 - j2)) == 1;
}

// const int N = consts::N_COLS * consts::N_ROWS;
constexpr int getN() {
    return consts::N_ROWS*consts::N_COLS;
}

Eigen::SparseMatrix<float> getA(int nRows, int nCols) {
    const int N = getN();
    Eigen::SparseMatrix<float> A(N, N);
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            A.coeffRef(i*nCols + j, i*nCols + j) = 4;

            if (j != 0) {
                A.coeffRef(i*nCols + j, i*nCols + j - 1) = -1;
            }
            if (j != nCols - 1) {
                A.coeffRef(i*nCols + j, i*nCols + j + 1) = -1;
            }
            if (i > 0) {
                A.coeffRef(i*nCols + j, (i-1)*nCols + j) = -1;
            }
            if (i < nRows - 1) {
                A.coeffRef(i*nCols + j, (i+1)*nCols + j) = -1;
            }
        }
    }
    A = A / (consts::delta_x*consts::delta_y);
    return A;
}

Eigen::VectorX<float> poisson(Eigen::SparseMatrix<float> A, Eigen::VectorX<float> f) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    solver.compute(A);
    return solver.solve(f);
    // Eigen::HouseholderQR<Eigen::MatrixX<float>> solver = A.householderQr();
    // Eigen::HouseholderQR<Eigen::MatrixX<float>> dec(A);
    // return solver.solve(f);
    // return A.householderQr().solve(f);
}

//  getF(Eigen::VectorX<float> u, int nRows, int nCols) {
//     for (int i = 0; i < nCols + 1; i++) {

//     }
// }
