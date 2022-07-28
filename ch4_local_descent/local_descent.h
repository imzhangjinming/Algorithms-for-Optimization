#ifndef __LOCAL_DESCENT_H__
#define __LOCAL_DESCENT_H__

#include<Eigen/Dense>

// using namespace Eigen;

// algorithm 4.2
// brief: backtracking line search
// param:
//      f: objective function
//      g: gradient of objective function
//      x: initial design point
//      d: descent direction
//      alpha: maximum step size
//      p: reduction factor, default = 0.5
//      beta: the first Wolfe condition parameter, default = 1e-4
// birth: created by ZJM on 20220728
// TODO: enable multivariate objective functions
template<typename T>
T backtrackingLineSearch(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>), Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>), Eigen::Matrix<T, Eigen::Dynamic, 1> x, Eigen::Matrix<T, Eigen::Dynamic, 1> d, T alpha, T p, T beta){
    T y = (*f)(x);
    Eigen::Matrix<T, Eigen::Dynamic, 1> dy = (*g)(x);
    while((*f)(x + alpha*d) > (y + beta * alpha * (dy.transpose() * d).value())){
        alpha *= p;
    }
    return alpha;
}

#endif