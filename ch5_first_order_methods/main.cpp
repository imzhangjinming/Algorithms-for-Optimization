#include<iostream>
#include<Eigen/Dense>
#include"first_order_method.h"
#include"local_descent.h"

using std::cout;
using std::endl;

// objective function
template<typename T>
T objective(Eigen::Matrix<T, Eigen::Dynamic, 1> x){
    // Rosenbrock function with a=1 and b=5
    // global minimum at [a, a^2]
    return (1 - x(0)) * (1 - x(0)) + 5 * (x(1) - x(0) * x(0)) * (x(1) - x(0) * x(0)); 
}

// gradient of objective function
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> gradient(Eigen::Matrix<T, Eigen::Dynamic, 1> x){
    return Eigen::Matrix<T, 2, 1>(20 * pow(x(0), 3) - 20 * x(0) * x(1) + 2 * x(0) - 2, 10 * (x(1)-x(0) * x(0)));
}



int main(void){
    // test algorithm 5.1
    cout << endl << "algorithm 5.1" << endl;
    Eigen::Matrix<double, 2, 1> x = {-1.0, -1.0};
    Eigen::Matrix<double, Eigen::Dynamic, 1> newX = gradientDescent<double>(objective<double>, gradient<double>, x, 1.25);
    cout << "newX = " << newX << endl;

    // test algorithm 5.2
    cout << endl << "algorithm 5.2" << endl;
    Eigen::Matrix<double, 2, 1> grad = {-44.0, -20.0};
    Eigen::Matrix<double, 2, 1> d = -grad;
    Eigen::Matrix<double, Eigen::Dynamic, 1> newX2 = conjugateGradient<double>(objective<double>, gradient<double>, x, d, grad);
    cout << "newX2 = " << newX2 << endl;

    // test algorithm 5.8
    cout << endl << "algorithm 5.2" << endl;
    Adam<double> adam;
    adam.alpha = 0.001;
    adam.gammaV = 0.9;
    adam.gammaS = 0.999;
    adam.epsilon = 1e-8;
    adam.k = 1;
    adam.v = Eigen::Matrix<double, 2, 1>::Zero();
    adam.s = Eigen::Matrix<double, 2, 1>::Zero();
    Eigen::Matrix<double, Eigen::Dynamic, 1> newX3 = adamDescent<double>(objective<double>, gradient<double>, x, adam);
    cout << "newX3 = " << newX3 << endl;

    return 0;
}