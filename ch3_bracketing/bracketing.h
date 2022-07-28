#ifndef __BRACKETING_H__
#define __BRACKETING_H__

#include<math.h>
#include<tuple>
#include<vector>
#include<limits>

// algorithm 3.1
// brief: find interval that contains minimum of unimodality function
// param: 
//      f: pointer of unimodality function
//      x: initial evaluation point
//      step: step size
//      k: expand factor
// birth: created by ZJM on 20220726
template<typename T>
std::tuple<T, T> bracketMinimum(T (*f)(T), T x, T step, T k){
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
            if(a < c){          
                return std::tuple<T, T>(a, c);
            }else{
                return std::tuple<T, T>(c, a);
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
//      f: pointer of unimodality function
//      a,b: interval contains minimum
//      n: max evaluation time of unimodality function
//      epsilon: 
// birth: created by ZJM on 20220727
template<typename T>
std::tuple<T, T> fibonacciSearch(T (*f)(T), T a, T b, int n, T epsilon){
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
    if(a < b){
        return std::tuple<T, T>(a, b);
    }
    else{
        return std::tuple<T, T>(b, a);
    }
}

// algorithm 3.3
// brief: golden section search
// param: 
//      f: pointer of unimodality function
//      a,b: interval contains minimum
//      n: max evaluation time of unimodality function
// birth: created by ZJM on 20220727
template<typename T>
std::tuple<T, T> goldenSectionSearch(T (*f)(T), T a, T b, int n){
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
    if(a < b){
        return std::tuple<T, T>(a, b);
    }
    else{
        return std::tuple<T, T>(b, a);
    } 
}

// algorithm 3.4
// brief: quadratic fit search
// param:
//      f: pointer of unimodality function
//      a,b,c: interval contains minimum, a < b < c
//      n: max evaluation time of unimodality function
// birth: created by ZJM on 20220727
template<typename T>
std::tuple<T, T, T> quadraticFitSearch(T (*f)(T), T a, T b, T c, int n){
    T ya = (*f)(a);
    T yb = (*f)(b);
    T yc = (*f)(c);
    T x, yx;
    for(int i = 0; i < n-3; i++){
        x = 0.5 * (ya * (b*b - c*c) + yb * (c*c - a*a) + yc * (a*a - b*b)) / (ya * (b-c) + yb * (c-a) + yc*(a-b));
        yx = (*f)(x);
        if(x > b){
            if(yx > yb){
                c = x;
                yc = yx;
            }
            else{
                a = b;
                ya = yb;
                b = x;
                yb = yx;
            }
        }
        else if(x < b){
            if(yx > yb){
                a = x;
                ya = yx;
            }
            else{
                c = b;
                yc = yb;
                b = x;
                yb = yx;
            }
        }        
    }
    return std::tuple<T, T, T>(a, b, c);
}

template<typename T>
struct Point2D{
    Point2D(){};

    template<typename T1>
    Point2D(T1 x, T1 y): x(x), y(y){}

    template<typename T2>
    Point2D<T2> &operator=(Point2D<T2> &rhs){ x=rhs.x; y=rhs.y; return *this;}

    T x;
    T y;
};

template<typename T>
Point2D<T> getSpIntersection(Point2D<T> A, Point2D<T> B, T l){
    T t = ((A.y - B.y) - l * (A.x - B.x)) / 2 / l;
    return Point2D<T>(A.x + t, A.y - t*l);
}

// algorithm 3.5
// brief: shubert-piyavskii method
// param: 
//      f: pointer of unimodality function
//      a,b: interval contains minimum
//      l:  Lipschitz constant
//      epsilon: termination tolerance
// birth: created by ZJM on 20220727
template<typename T>
Point2D<T> shubertPiyavskii(T (*f)(T), T a, T b, T l, T epsilon){
    T m = (a + b) / 2;
    Point2D<T> A(a, (*f)(a));
    Point2D<T> M(m, (*f)(m));
    Point2D<T> B(b, (*f)(b));
    Point2D<T> AM = getSpIntersection(A, M, l);
    Point2D<T> MB = getSpIntersection(M, B, l);
    std::vector<Point2D<T>> pts{A, AM, M, MB, B};
    T diff = std::numeric_limits<T>::infinity();
    Point2D<T> p, pPrev, pNext;
    int i;
    while(diff > epsilon){
        i = 0;
        for(int idx = 0; idx < pts.size(); idx++){
            if(pts[idx].y < pts[i].y) i = idx;
        }
        p.x = pts[i].x;
        p.y = (*f)(p.x);
        diff = p.y - pts[i].y;

        pPrev = getSpIntersection(pts[i-1], p, l);
        pNext = getSpIntersection(p, pts[i+1], l);

        pts[i] = p;
        pts.insert(pts.begin() + i + 1, pNext);
        pts.insert(pts.begin() + i, pPrev);
    }
    i = 0;
    for(int idx = 0; idx < pts.size(); idx+=2){
        if(pts[idx].y < pts[i].y) i = idx;
    }
    return pts[i];
}

// algorithm 3.6
// brief: bisection algorithm
// param:
//      g: derivative of univariate function
//      a,b: initial interval
//      epsilon: interval width tolerance
// birth: created by ZJM on 20220728
template<typename T>
std::tuple<T, T> bisection(T (*g)(T), T a, T b, T epsilon){
    if(a > b){ a = a + b; b = a - b; a = a - b;}
    T ya = (*g)(a);
    T yb = (*g)(b);
    if(ya == 0) b = a;
    if(yb == 0) a = b;
    T m, ym;
    while(b - a > epsilon){
        m = (a + b) / 2;
        ym = (*g)(m);
        if(ym == 0){
            a = m;
            b = m;
        }
        else if(ym * ya > 0){
            a = m;
        }
        else if(ym * yb > 0){
            b = m;
        }
    }
    return std::tuple<T, T>(a, b);
}

#endif