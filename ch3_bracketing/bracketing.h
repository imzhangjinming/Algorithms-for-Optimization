#ifndef __BRACKETING_H__
#define __BRACKETING_H__

#include<math.h>

template<typename T>
struct Interval
{
    T lowerBound;
    T upperBound;
};

// algorithm 3.1
// brief: 计算单峰函数极小值点所在的区间
// param: 
//      f: 指向单峰函数的指针
//      x: 初始点
//      step: 寻找极小值点时采用的步长
//      k: 步长增长的倍数
// birth: created by ZJM on 20220726
template<typename T>
Interval<T> bracketMinimum(T (*f)(T), T x, T step, T k){
    T a = x;
    T ya = (*f)(a);
    T b = x + step;
    T yb = (*f)(b);
    if (yb > ya){
        step = -step;
        b = x + step;
        yb = (*f)(b);
    }
    T c, yc;
    while(1){
        c = b + step;
        yc = (*f)(c);
        if(yc > yb){
            Interval<T> interval;
            if(a < c){          
                interval.lowerBound = a;
                interval.upperBound = c;
                return interval;
            }else{
                interval.lowerBound = c;
                interval.upperBound = a;
                return interval;
            }
        }
        a = b;
        ya = yb;
        b = c;
        yb = yc;
        step *= k;
    }
}

// algorithm 3.2
// brief: fibonacci search
// param: 
//      f: 指向单峰函数的指针
//      a,b: 包含函数极小值点的区间
//      n: 调用单峰函数的次数
//      epsilon: 
// birth: created by ZJM on 20220727
template<typename T>
Interval<T> fibonacciSearch(T (*f)(T), T a, T b, int n, T epsilon){
    T s = (1 - sqrt(5)) / (1 + sqrt(5));
    T phi = (1 + sqrt(5)) / 2;
    T rho = 1 / (phi * (1 - pow(s, n+1)) / (1 - pow(s, n)));
    T d = rho * b + (1 - rho) * a;
    T yd = (*f)(d);
    T c, yc;
    for(int i = 0; i < n-1; i++){
        if(i == n-2){
            c = epsilon * a + (1 - epsilon) * d;
        }
        else{
            c = rho * a + (1 - rho) * b;
        }
        yc = (*f)(c);
        if(yc < yd){
            b = d;
            d = c;
            yd = yc;
        }
        else{
            a = b;
            b = c;
        }
        rho = 1 / (phi * (1 - pow(s, n-i)) / (1 - pow(s, n-i-1)));
    }
    Interval<T> interval;
    if(a < b){
        interval.lowerBound = a;
        interval.upperBound = b;
    }
    else{
        interval.lowerBound = b;
        interval.upperBound = a;
    }
    return interval;
}

// algorithm 3.3
// brief: golden section search
// param: 
//      f: 指向单峰函数的指针
//      a,b: 包含函数极小值点的区间
//      n: 调用单峰函数的次数
// birth: created by ZJM on 20220727
template<typename T>
Interval<T> goldenSectionSearch(T (*f)(T), T a, T b, int n){
    T phi = (1 + sqrt(5)) / 2;
    T rho = phi - 1;
    T d = rho * b + (1 - rho) * a;
    T yd = (*f)(d);
    T c, yc;
    for(int i = 0; i < n-1; i++){
        c = rho * a + (1 - rho) * b;
        yc = (*f)(c);
        if(yc < yd){
            b = d;
            d = c;
            yd = yc;
        }
        else{
            a = b;
            b = c;
        }
    }
    Interval<T> interval;
    if(a < b){
        interval.lowerBound = a;
        interval.upperBound = b;
    }
    else{
        interval.lowerBound = b;
        interval.upperBound = a;
    }
    return interval;    
}


#endif