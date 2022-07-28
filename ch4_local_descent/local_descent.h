#ifndef __LOCAL_DESCENT_H__
#define __LOCAL_DESCENT_H__

#include<Eigen/Dense>
#include<limits>
#include<cmath>
#include<cstdlib>

// algorithm 4.2
// brief: backtracking line search
// param:
//      f: objective function
//      g: gradient of objective function
//      x: initial design point
//      d: descent direction
//      alpha: maximum step size
//      p: reduction factor
//      beta: the first Wolfe condition parameter
// return: step size
// birth: created by ZJM on 20220728
template<typename T>
T backtrackingLineSearch(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>), Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>), Eigen::Matrix<T, Eigen::Dynamic, 1> x, Eigen::Matrix<T, Eigen::Dynamic, 1> d, T alpha, T p, T beta){
    T y = (*f)(x);
    Eigen::Matrix<T, Eigen::Dynamic, 1> dy = (*g)(x);
    while((*f)(x + alpha*d) > (y + beta * alpha * (dy.transpose() * d).value())){
        alpha *= p;
    }
    return alpha;
}

// algorithm 4.3
// brief: strong backtracking approximate line search for satisfying the strong Wolfe conditions
// param:
//      f: objective function
//      g: gradient of objective function
//      x: initial design point
//      d: descent direction
//      alpha: initial step size
//      beta: Wolfe condition paramter
//      sigma: Wolfe condition paramter
// return: step size
// birth: created by ZJM on 20220728
template<typename T>
T strongBacktracking(T (*f)(Eigen::Matrix<T, Eigen::Dynamic, 1>), \
                        Eigen::Matrix<T, Eigen::Dynamic, 1> (*g)(Eigen::Matrix<T, Eigen::Dynamic, 1>), \
                        Eigen::Matrix<T, Eigen::Dynamic, 1> x, \
                        Eigen::Matrix<T, Eigen::Dynamic, 1> d, \
                        T alpha, T beta, T sigma){
    T y0 = (*f)(x), dy0 = ((*g)(x).transpose() * d).value();
    T yPrev = std::nan("1"), alphaPrev = 0;
    T alphaLow = std::nan("1"), alphaHigh = std::nan("1");

    // bracket phase
    T y, dy;
    while(1){
        y = (*f)(x + alpha * d);
        if((y > y0 + beta * alpha * dy0) || (!std::isnan(yPrev)&&(y >= yPrev))){
            alphaLow = alphaPrev;
            alphaHigh = alpha;
            break;
        }
        dy = ((*g)(x + alpha * d).transpose() * d).value();
        if(std::abs(dy) <= -sigma*dy0){
            return alpha;
        }
        else if(dy >= 0){
            alphaLow = alpha;
            alphaHigh = alphaPrev;
            break;
        }
        yPrev = y;
        alphaPrev = alpha;
        alpha *= 2;
    }

    // zoom phase
    T yLow = (*f)(x + alphaLow * d);
    while(1){
        alpha = (alphaLow + alphaHigh) / 2;
        y = (*f)(x + alpha * d);
        if(y > y0 + beta*alpha*dy0 || y >= yLow){
            alphaHigh = alpha;
        }
        else{
            dy = ((*g)(x + alpha*d).transpose()*d).value();
            if(std::abs(dy) <= -sigma*dy0){
                return alpha;
            }
            else if(dy*(alphaHigh - alphaLow) >= 0){
                alphaHigh = alphaLow;
            }
            alphaLow = alpha;
        }
    }
}


#endif