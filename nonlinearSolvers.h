//
// Created by Kuboss on 22.05.2023.
//

#ifndef MO10_NONLINEARSOLVERS_H
#define MO10_NONLINEARSOLVERS_H
#include <iostream>
#include <math.h>
#include <functional>
using namespace std;
#define ERR 1.3e-6
#define FUNERR 1.3e-6


double picard(function<double(double)>fx, double initialGuess, int MaxIteracji);
double bisekcja(function<double(double)>fx, double a, double b, int MaxIteracji);
double newton(function<double(double)>fx, function<double(double)>derivative, double initialGuess, int MaxIteracji);
double sieczne(function<double(double)>fx, double initialGuess1, double initialGuess2, int MaxIteracji);
#endif //MO10_NONLINEARSOLVERS_H
