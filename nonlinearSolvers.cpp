//
// Created by Kuboss on 22.05.2023.
//
#include "nonlinearSolvers.h"
double picard(function<double(double)>fx, double initialGuess, int MaxIteracji){
    double xn=initialGuess, xn1;
    for(int i=0; i<MaxIteracji; i++){
        xn1=fx(xn)+xn;
        if(i+1==MaxIteracji){
            //cout << "Zakończono z powodu wyczerpania liczby iteracji";
            return numeric_limits<double>::quiet_NaN();
        }
        //cout << "Iteracja nr " << i+1 << " xn+1= " << xn1 << ". Estymator bledu: " << xn1-xn << ". Residuum rownania: "<< fx(xn1) <<  endl;
        if(abs(xn1-xn)<=ERR && abs(fx(xn1))<=FUNERR){
            //cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn" << endl;
            return xn1;
        }
        xn=xn1;
    }
    //cout << "Pierwiastek: " << xn1 << endl;
}

double bisekcja(function<double(double)>fx, double a, double b, int MaxIteracji){
    double x;
    if(fx(a)<0 && fx(b)<0 || fx(a)>0 && fx(b)>0){
        cout << "Brak pierwiastków w podanym przedziale"<< endl;
    }
    if(a>b){
        swap(a,b);
    }
    for(int i=0; i<MaxIteracji; i++){
        x=(a+b)/2.0;
        if(i+1==MaxIteracji){
            cout << "Zakonczono z powodu wyczerpania liczby iteracji"<< endl;
            return numeric_limits<double>::quiet_NaN();
        }
        //cout << "Iteracja nr " << i+1 << "  Przedzial: ["<< a << "," << b << "]" << " xn= " << x  << ". Estymator bledu: " << (b-a)/2 << ". Residuum rownania: "<< fx(x) << endl;
        if(abs((b-a)/2)<=ERR && abs(fx(x))<=FUNERR){
            //cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn" << endl;
            return x;
        }


        if(x<0 && b>0 || x>0 && b<0){
            a=x;
        }else{
            b=x;
        }
    }
    //cout << "Pierwiastek: " << x << endl;
}

double newton(function<double(double)>fx, function<double(double)>derivative, double initialGuess, int MaxIteracji){
    double xn=initialGuess, xn1;
    for(int i=0; i<MaxIteracji; i++){
        xn1=xn - fx(xn)/derivative(xn);
        if(i+1==MaxIteracji){
            cout << "Zakonczono z powodu wyczerpania liczby iteracji"<< endl;
            return numeric_limits<double>::quiet_NaN();
        }
        //cout << "Iteracja nr " << i+1 << " xn+1= " << xn1 << ". Estymator bledu: " << xn1-xn << ". Residuum rownania: "<< fx(xn1) <<  endl;
        if(abs(xn1-xn)<=ERR && abs(fx(xn1))<=FUNERR){
            //cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn" << endl;
            return xn1;
        }
        xn=xn1;
    }
    //cout << "Pierwiastek: " << xn1 << endl;
}

double sieczne(function<double(double)>fx, double initialGuess1, double initialGuess2, int MaxIteracji){
    double xn=initialGuess1, xn1 = initialGuess2, xn2;
    for(int i=0; i<MaxIteracji; i++){
        xn2=xn1 - fx(xn1)/((fx(xn1)-fx(xn))/(xn1-xn));
        if(i+1==MaxIteracji){
            //cout << "Zakonczono z powodu wyczerpania liczby iteracji"<< endl;
            return numeric_limits<double>::quiet_NaN();
        }
        //cout << "Iteracja nr " << i+1 << " xn+1= " << xn2 << ". Estymator bledu: " << xn2-xn1 << ". Residuum rownania: "<< fx(xn2) <<  endl;
        if(abs(xn2-xn1)<=ERR && abs(fx(xn2))<=FUNERR){
            //cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn" << endl;
            return xn2;
        }
        xn=xn1;
        xn1=xn2;
    }
    //cout << "Pierwiastek: " << xn2 << endl;
}

