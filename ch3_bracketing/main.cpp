#include"bracketing.h"
#include<iostream>

using std::cout;
using std::endl;

// Unimodality function
template<typename T>
T unimodalityFunc(T x){
    return x * x - 2 * x + 1;
}

int main(void){
    // test algorithm 3.1
    cout << "algorithm 3.1" << endl;
    Interval<double> interval1 = bracketMinimum<double>(unimodalityFunc<double>, 0.0, 1e-2, 2.0);
    cout << "lower bound: " << interval1.lowerBound << endl;
    cout << "upper bound: " << interval1.upperBound << endl;

    Interval<float> interval2 = bracketMinimum<float>(unimodalityFunc<float>, 0.0, 1e-2, 2.0);
    cout << "lower bound: " << interval1.lowerBound << endl;
    cout << "upper bound: " << interval1.upperBound << endl;

    Interval<int> interval3 = bracketMinimum<int>(unimodalityFunc<int>, 0, 1, 2);
    cout << "lower bound: " << interval3.lowerBound << endl;
    cout << "upper bound: " << interval3.upperBound << endl;
    
    // test algorithm 3.2
    cout << endl << "algorithm 3.2" << endl;
    Interval<double> interval4 = fibonacciSearch<double>(unimodalityFunc<double>, interval1.lowerBound, interval1.upperBound, 8, 0.01l);
    cout << "lower bound: " << interval4.lowerBound << endl;
    cout << "upper bound: " << interval4.upperBound << endl;

    Interval<float> interval5 = fibonacciSearch<float>(unimodalityFunc<float>, interval2.lowerBound, interval2.upperBound, 8, 0.01f);
    cout << "lower bound: " << interval5.lowerBound << endl;
    cout << "upper bound: " << interval5.upperBound << endl;       

    // test algorithm 3.3
    cout << endl << "algorithm 3.3" << endl;
    Interval<double> interval6 = goldenSectionSearch<double>(unimodalityFunc<double>, interval1.lowerBound, interval1.upperBound, 8);
    cout << "lower bound: " << interval6.lowerBound << endl;
    cout << "upper bound: " << interval6.upperBound << endl;

    Interval<float> interval7 = fibonacciSearch<float>(unimodalityFunc<float>, interval2.lowerBound, interval2.upperBound, 8, 0.01f);
    cout << "lower bound: " << interval7.lowerBound << endl;
    cout << "upper bound: " << interval7.upperBound << endl;  

    

    return 0;
}