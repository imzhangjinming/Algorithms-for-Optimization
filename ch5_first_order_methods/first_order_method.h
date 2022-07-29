#ifndef __FIRST_ORDER_METHOD_H__
#define __FIRST_ORDER_METHOD_H__

#include<Eigen/Dense>
#include<cmath>
#include"local_descent.h"

// algorithm 5.1
// brief: gradient descent method
// param:
//      f: objective function
//      g: gradient of objective function
//      x: present design point
//      alpha: step size(learning rate)
// return: new design point
// birth: created by ZJM on 20220729
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> gradientDescent(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>),
                                            Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>),
                                            Eigen::Matrix<T, Eigen::Dynamic, 1> x,
                                            T alpha){
    return x - alpha * (*g)(x) / (*g)(x).norm();
}

// algorithm 5.2
// brief: conjugate gradient descent
// param:
//      f: objective function
//      g: gradient of objective function
//      x: present design point
//      d: previous search direction
//      grad: previous gradient
// return: new design point
// birth: created by ZJM on 20220729
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> conjugateGradient(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> x,\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> d,\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> grad){
    Eigen::Matrix<T, Eigen::Dynamic, 1> gradPresent = (*g)(x);
    T beta = (gradPresent.transpose() * (gradPresent - grad)).value() / grad.squaredNorm();
    beta = beta > 0 ? beta : 0;
    Eigen::Matrix<T, Eigen::Dynamic, 1> dPresent = -gradPresent + beta * d;
    T alpha = backtrackingLineSearch<T>(f, g, x, dPresent, 10.0, 0.5, 1e-4);
    Eigen::Matrix<T, Eigen::Dynamic, 1> xNew = x + alpha * dPresent;
    return xNew;
}

// algorithm 5.3 
// brief: the momentum method for accelerated descent
// param:
//      f: objective function
//      g: gradient of objective function
//      x: present design point
//      alpha: learning rate
//      beta: momentum decay
//      v: momentum
// return: new design point
// birth: created by ZJM on 20220729
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> momentumDescent(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> x,\
                                                T alpha,\
                                                T beta,\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> &v){
    Eigen::Matrix<T, Eigen::Dynamic, 1> grad = (*g)(x);
    v = beta * v - alpha * grad;
    return x + v;
}

// algorithm 5.4
// brief: Nesterov’s momentum method of accelerated descent
// param:
//      f: objective function
//      g: gradient of objective function
//      x: present design point
//      alpha: learning rate
//      beta: momentum decay
//      v: momentum
// return: new design point
// birth: created by ZJM on 20220729
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> nesterovMomentum(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> x,\
                                                T alpha,\
                                                T beta,\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> &v){
    v = beta * v - alpha * (*g)(x + beta * v);
    return x + v;
}

template<typename T>
struct Adam{
    T alpha; // learning rate
    T gammaV; // decay
    T gammaS; // decay
    T epsilon; // small value
    int k; // step counter
    Eigen::Matrix<T, Eigen::Dynamic, 1> v; // 1st moment estimate
    Eigen::Matrix<T, Eigen::Dynamic, 1> s; // 2nd moment estimate
};

// algorithm 5.8
// brief: the Adam accelerated descent method
// param:
//      f: objective function
//      g: gradient of objective function
//      x: present design point
//      adam: instance of struct Adam
// return: new design point
// birth: created by ZJM on 20220729
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> adamDescent(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>),\
                                                Eigen::Matrix<T, Eigen::Dynamic, 1> x,\
                                                Adam<T> &adam){
    T alpha = adam.alpha, gammaV = adam.gammaV, gammaS = adam.gammaS, epsilon = adam.epsilon;
    int k = adam.k;
    Eigen::Matrix<T, Eigen::Dynamic, 1> grad = (*g)(x);
    adam.v = gammaV * adam.v + (1 - gammaV) * grad;
    adam.s = gammaS * adam.s + (1 - gammaS) * grad.cwiseProduct(grad);
    adam.k += 1;
    Eigen::Matrix<T, Eigen::Dynamic, 1> vHat = adam.v / (1 - pow(gammaV, k));
    Eigen::Matrix<T, Eigen::Dynamic, 1> sHat = adam.s / (1 - pow(gammaS, k));
    return x - alpha * vHat.cwiseQuotient(sHat.cwiseSqrt()); //+ epsilon
}

#endif