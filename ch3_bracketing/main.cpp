#include"bracketing.h"
#include<iostream>

using std::cout;
using std::endl;
using std::tuple;
using std::get;

// Unimodality function
template<typename T>
T unimodalityFunc(T x){
    return exp(-x + 1) + x - 1;  //;x * x - 2 * x + 1
}

int main(void){
    // test algorithm 3.1
    cout << "algorithm 3.1" << endl;
    tuple<double, double> interval1 = bracketMinimum<double>(unimodalityFunc<double>, 0.0, 1e-2, 2.0);
    cout << "lower bound: " << get<0>(interval1) << endl;
    cout << "upper bound: " << get<1>(interval1) << endl;

    tuple<float, float> interval2 = bracketMinimum<float>(unimodalityFunc<float>, 0.0, 1e-2, 2.0);
    cout << "lower bound: " << get<0>(interval2) << endl;
    cout << "upper bound: " << get<1>(interval2) << endl;

    tuple<int, int> interval3 = bracketMinimum<int>(unimodalityFunc<int>, 0, 1, 2);
    cout << "lower bound: " << get<0>(interval3) << endl;
    cout << "upper bound: " << get<1>(interval3) << endl;
    
    // test algorithm 3.2
    cout << endl << "algorithm 3.2" << endl;
    tuple<double, double> interval4 = fibonacciSearch<double>(unimodalityFunc<double>, get<0>(interval1), get<1>(interval1), 8, 0.01);
    cout << "lower bound: " << get<0>(interval4) << endl;
    cout << "upper bound: " << get<1>(interval4) << endl;

    tuple<float, float> interval5 = fibonacciSearch<float>(unimodalityFunc<float>, get<0>(interval2), get<1>(interval2), 8, 0.01f);
    cout << "lower bound: " << get<0>(interval5) << endl;
    cout << "upper bound: " << get<1>(interval5) << endl;       

    // test algorithm 3.3
    cout << endl << "algorithm 3.3" << endl;
    tuple<double, double> interval6 = goldenSectionSearch<double>(unimodalityFunc<double>, get<0>(interval1), get<1>(interval1), 8);
    cout << "lower bound: " << get<0>(interval6) << endl;
    cout << "upper bound: " << get<1>(interval6) << endl;

    tuple<float, float> interval7 = goldenSectionSearch<float>(unimodalityFunc<float>, get<0>(interval2), get<1>(interval2), 8);
    cout << "lower bound: " << get<0>(interval7) << endl;
    cout << "upper bound: " << get<1>(interval7) << endl;  

    // test algorithm 3.4
    cout << endl << "algorithm 3.4" << endl;
    tuple<double, double, double> interval8 = quadraticFitSearch<double>(unimodalityFunc<double>, get<0>(interval1), (get<0>(interval1) + get<1>(interval1)) / 2, get<1>(interval1), 8);
    double a1 = get<0>(interval8);
    double b1 = get<1>(interval8);
    double c1 = get<2>(interval8);
    double ya1 = unimodalityFunc<double>(a1);
    double yb1 = unimodalityFunc<double>(b1);
    double yc1 = unimodalityFunc<double>(c1);
    double x_star1 = 0.5 * (ya1 * (b1*b1 - c1*c1) + yb1 * (c1*c1 - a1*a1) + yc1 * (a1*a1 - b1*b1)) / (ya1 * (b1-c1) + yb1 * (c1-a1) + yc1 * (a1-b1));
    cout << "x_star: " << x_star1 << endl;

    tuple<float, float, float> interval9 = quadraticFitSearch<float>(unimodalityFunc<float>, get<0>(interval2), (get<0>(interval2) + get<1>(interval2)) / 2, get<1>(interval2), 8);
    float a = get<0>(interval9);
    float b = get<1>(interval9);
    float c = get<2>(interval9);
    float ya = unimodalityFunc<float>(a);
    float yb = unimodalityFunc<float>(b);
    float yc = unimodalityFunc<float>(c);
    float x_star = 0.5 * (ya * (b*b - c*c) + yb * (c*c - a*a) + yc * (a*a - b*b)) / (ya * (b-c) + yb * (c-a) + yc * (a-b));
    cout << "x_star: " << x_star << endl;

    // test algorithm 3.5
    cout << endl << "algorithm 3.5" << endl;
    Point2D<double> p1 = shubertPiyavskii(unimodalityFunc<double>, get<0>(interval1), get<1>(interval1), 5.0, 0.1);
    cout << "x: " << p1.x << ", " << "y: " << p1.y << endl;

    Point2D<float> p2 = shubertPiyavskii<float>(unimodalityFunc<float>, get<0>(interval2), get<1>(interval2), 5.0, 0.1);
    cout << "x: " << p2.x << ", " << "y: " << p2.y << endl; 

    return 0;
}