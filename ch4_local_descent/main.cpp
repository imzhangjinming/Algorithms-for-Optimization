#include<iostream>
#include<Eigen/Dense>
#include"local_descent.h"

using std::cout;
using std::endl;

// multivariate objective function
template<typename T>
T objective(Eigen::Matrix<T, Eigen::Dynamic, 1> x){
    return x(0) * x(0) + x(0) * x(1) + x(1) * x(1);
}

// derivate of objective function
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> derivative(Eigen::Matrix<T, Eigen::Dynamic, 1> x){
    return Eigen::Matrix<T, 2, 1>(2 * x(0) + x(1), x(0) + 2 * x(1));
}


int main(void){
    // test algorithm 4.2
    cout << endl << "algorithm 4.2" << endl;
    Eigen::Matrix<double, 2, 1> designPoint = {1, 2};
    Eigen::Matrix<double, 2, 1> direction = {-1, -1};

    double alpha = backtrackingLineSearch<double>(objective<double>, derivative<double>, \
                                                    designPoint, direction, 10.0, 0.5, 1e-4);
    cout << "alpha = " << alpha << endl;

    return 0;
}