#include <iostream>
#include <functional>
#include <fstream>
#include <math.h>
#include "nonlinearSolvers.h"


using namespace std;

double analityczne(double t){
    return (1. - exp(-10.*(t + atan(t))));
}



void BME(function<double(double, double)>f, double dt, double yk, ofstream &plik){
    double yk1;
    plik << 0 << " " << yk << '\n';
    double t = 0+dt;
    while(true){
        yk1 = yk + f(t, yk)*dt;
        plik << t << " " << yk1 << '\n';
        if(abs(yk-yk1)<10e-6 || abs(yk-yk1)>10e3){
            break;
        }
        yk = yk1;
        t+=dt;
    }
}

void PME(double dt, double yk, ofstream &plik){
    double yk1;
    plik << 0 << " " << yk << '\n';
    double t = 0+dt;
    while(true){
        yk1 = sieczne([yk,dt, t](double x)->double{return (x-yk)/dt+((10.0*(t+dt)*(t+dt)+20.0)/(t*t+1))*(x-1);}, yk, yk+dt, 100);
        plik << t << " " << yk1 << '\n';
        if(abs(yk-yk1)<10e-6){
            break;
        }
        yk = yk1;
        t+=dt;
    }
}

void PMT(double dt, double yk, ofstream &plik){
    double yk1;
    plik << 0 << " " << yk << '\n';
    double t = 0+dt;
    while(true){
        yk1 = sieczne([yk,dt, t](double x)->double{return (x-yk)/dt+(((10.0*(t+dt)*(t+dt)+20.0)/(t*t+1))*(x-1) + (10.0*t*t+20.0)/(t*t+1)*(yk-1))/2;}, yk, yk+dt, 100);
        plik << t << " " << yk1 << '\n';
        if(abs(yk-yk1)<10e-6){
            break;
        }
        yk = yk1;
        t+=dt;
    }
}

void BME_error(function<double(double, double)>f, double dt, double yk, ofstream &plik){
    double yk1, blad, maxBlad=0;
    //plik << 0 << " " << yk << '\n';
    double t = 0+dt;
    while(true){
        yk1 = yk + f(t, yk)*dt;

        if(abs(yk-yk1)<10e-6){
            break;
        }
        blad = fabs(yk1- analityczne(t+dt));
        if(blad>maxBlad)
            maxBlad = blad;
        yk = yk1;
        t+=dt;
    }
    plik << log10(dt) << " " << log10(maxBlad)<< '\n';
}

void PME_error(double dt, double yk, ofstream &plik){
    double yk1, blad, maxBlad=0;
    //plik << 0 << " " << yk << '\n';
    double t = 0+dt;
    while(true){
        yk1 = sieczne([yk,dt, t](double x)->double{return (x-yk)/dt+((10.0*(t+dt)*(t+dt)+20.0)/(t*t+1))*(x-1);}, yk, yk+dt, 100);
        //plik << t << " " << yk1 << '\n';
        if(abs(yk-yk1)<10e-6){
            break;
        }
        blad = fabs(yk1- analityczne(t+dt));
        if(blad>maxBlad)
            maxBlad = blad;
        yk = yk1;
        t+=dt;
    }
    plik << log10(dt) << " " << log10(maxBlad)<< '\n';
}

void PMT_error(double dt, double yk, ofstream &plik){
    double yk1, blad, maxBlad=0;
    //plik << 0 << " " << yk << '\n';
    double t = 0+dt;
    while(true){
        yk1 = sieczne([yk,dt, t](double x)->double{return (x-yk)/dt+(((10.0*(t+dt)*(t+dt)+20.0)/(t*t+1))*(x-1) + (10.0*t*t+20.0)/(t*t+1)*(yk-1))/2;}, yk, yk+dt, 100);
        //plik << t << " " << yk1 << '\n';
        if(abs(yk-yk1)<10e-6){
            break;
        }
        blad = fabs(yk1- analityczne(t));
        if(blad>maxBlad)
            maxBlad = blad;
        yk = yk1;
        t+=dt;
    }
    plik << log10(dt) << " " << log10(maxBlad)<< '\n';
}

int main() {
    double dt = 0.001;
    ofstream plik("wynikBME.txt");
    BME([](double t, double y)->double{return -((10.0*t*t+20.0)/(t*t+1))*(y-1);}, dt, 0, plik);
    plik.close();

    plik.open("wynikBMEniestabilny.txt");
    BME([](double t, double y)->double{return -((10.0*t*t+20.0)/(t*t+1))*(y-1);}, 0.15, 0, plik);
    plik.close();

    plik.open("wynikPME.txt");
    PME(dt, 0, plik);
    plik.close();

    plik.open("wynikPMT.txt");
    PMT(dt, 0, plik);
    plik.close();


    plik.open("wynikAnalityczny.txt");
    double t = 0+dt;
    do{
        plik << t << " " << analityczne(t) << '\n';
        t+=dt;
    }while(abs(analityczne(t)-analityczne(t-dt))>10e-6);
    plik.close();

    plik.open("wynikAnalitycznydoNiestabilnego.txt");
    t = 0+0.01;
    do{
        plik << t << " " << analityczne(t) << '\n';
        t+=0.01;
    }while(abs(analityczne(t)-analityczne(t-dt))>10e-25);
    plik.close();

    plik.open("BMEBlad.txt");
    for(double i=0.1;i>10e-15;i/=2){
        BME_error([](double t, double y)->double{return -((10.0*t*t+20.0)/(t*t+1))*(y-1);}, i, 0, plik);
    }
    plik.close();

    plik.open("PMEBlad.txt");
    for(double i=0.1;i>10e-15;i/=2){
        PME_error( i, 0, plik);
    }
    plik.close();

    plik.open("PMTBlad.txt");
    for(double i=0.1;i>10e-15;i/=2){
        PMT_error( i, 0, plik);
    }
    plik.close();

    return 0;
}
