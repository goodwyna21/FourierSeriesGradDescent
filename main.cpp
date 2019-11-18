#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

struct Descender{
    const unsigned int N; //Num Terms
    const unsigned int D; //Num Data Points
    double* d; //Data Points x values
    double* y; //Data Points y values

    Descender(const unsigned int terms, const unsigned int points, double* datax, double* datay)
    :N(terms),D(points),d(datax),y(datay){}
    Descender():N(0),D(0){}
    ~Descender(){}

    static void genPoints(uint length, double x0, double x1, double (*f)(double), double* d, double* y){
        for(int i = 0; i < length; i++){
            d[i] = x0 + ((1.0*i/(length-1)) * (x1 - x0));
            y[i] = (*f)(d[i]);
        }
    }

    double k(double x){} //evauluate
    double g(){} //cost
    void iterate(const double stepsize){}
};

struct Polynomial : public Descender{
    double* c; //coefficients

    Polynomial(unsigned int degree, unsigned int numdata,double* datain, double* dataout, double minc, double maxc)
    :Descender(degree,numdata,datain,dataout){
        c = new double[N+1];
        for(int n = 0; n <= N; n++){
            c[n] = ((double)rand()/RAND_MAX)*(maxc-minc) + minc;
        }
    }
    Polynomial(){}
    ~Polynomial(){}

    double k(double x){
        double sum = 0;
        for(int n = 0; n <= N; n++){
            sum += c[n] * pow(x,n);
        }
        return sum;
    }

    double g(){
        double sum = 0;
        for(int i = 0; i < D; i++){
            sum += (1.0/D)*pow(k(d[i]) - y[i], 2);
        }
        return sum;
    }

    void iterate(const double stepsize){
        double* newc = new double[N];
        for(int n = 0; n <= N; n++){
            newc[n] = c[n] - stepsize*dgdcn(n);
        }
        for(int n = 0; n <= N; n++){
            c[n] = newc[n];
        }
    }



    double dgdcn(unsigned int n){
        double sum = 0;
        for(int i = 0; i < D; i++){
            sum += (2.0/D) * (k(d[i]) - y[i]) * pow(d[i],n);
        }
        return sum;
    }
};

struct Fourier : public Descender{
    double* P; //Periods
    double* a; //Coefficients

    Fourier(unsigned int size, unsigned int numdata,double* datain, double* dataout, double mina, double maxa, double p0, double p1)
    :Descender(size,numdata,datain,dataout){
        a = new double[N];
        P = new double[N];
        for(int n = 0; n < N; n++){
            a[n] = ((double)rand()/RAND_MAX)*(maxa-mina) + mina;
            P[n] = p0 + ((1.0*n/(N-1)) * (p1 - p0));
        }
    }
    Fourier(){}
    ~Fourier(){}

    double k(double x){
        double sum = 0;
        for(int n = 0; n < N; n++){
            sum += a[n] * sin(x*P[n]);
        }
        return sum;
    }

    double g(){
        double sum = 0;
        for(int i = 0; i < D; i++){
            sum += (1.0/D)*pow(k(d[i]) - y[i], 2);
        }
        return sum;
    }

    void iterate(const double stepsize){
        double* newa = new double[N];
        for(int n = 0; n < N; n++){
            newa[n] = a[n] - stepsize*dgdan(n);
        }
        for(int n = 0; n < N; n++){
            a[n] = newa[n];
        }
        delete[] newa;
    }



    double dgdan(unsigned int n){
        double sum = 0;
        for(int i = 0; i < D; i++){
            sum += (2.0/D) * (k(d[i]) - y[i]) * sin(d[i] * P[n]);
        }
        return sum;
    }
};

// -- EXAMPLE FUNCTIONS --
double squarewave(double x){
    return 2*((int)x%2) - 1;
}

double smoothfunct(double x){
    return 0.5*sin(x) + 0.5*sin(x/2) + sin(-2/3.0*x);
}

double polynomialfunct(double x){
    return -0.8*pow(x,3) - 2*pow(x,2) - x + 0.6;
}

void PolynomialDescent(unsigned int size,unsigned int numpoints,unsigned int iterations,double tolerance,double stepsize,double minc,double maxc,double(*f)(double),double x0,double x1){
    double xs[numpoints];
    double ys[numpoints];
    Descender::genPoints(numpoints,x0,x1,f,xs,ys);

    Polynomial d(size,numpoints,xs,ys,minc,maxc);

    cout << "Polynomial Gradient Decent" <<
    "\nDegree     : " << size <<
    "\nNum Points : " << numpoints <<
    "\nIterations : " << iterations <<
    "\nTolerance  : " << tolerance <<
    "\nStep Size  : " << stepsize <<
    "\nmin c      : " << minc <<
    "\nmax c      : " << maxc <<
    "\nx0         : " << x0 <<
    "\nx1         : " << x1 << "\n\n";

    cout << "Initial coefficients:\n";
    for(int n = 0; n < size+1; n++){
        cout << d.c[n] << ", ";
    }
    cout << "\n\nData x's:\n";
    for(int i = 0; i < numpoints; i++){
        cout << xs[i] << ", ";
    }
    cout << "\n\nData y's:\n";
    for(int i = 0; i < numpoints; i++){
        cout << ys[i] << ", ";
    }
    cout << "\n\nIterations:\n";

    double cost;
    for(int i = 0; i < iterations; i++){
        cost = d.g();
        if(i%10 == 0){
            cout << i << ": " << cost << "\n";
        }
        d.iterate(stepsize);
        if(cost <= tolerance){
            cout << "\nConverged at epoch " << i << "\n";
            break;
        }
    }
    cout << "Final: " << d.g() << "\n";

    cout << "\n\ncoefficients:\n";
    for(int n = 0; n < size+1; n++){
        cout << d.c[n] << ", ";
    }
    cout << "\n";
}

void FourierDescent(unsigned int size,unsigned int numpoints,unsigned int iterations,double tolerance,double stepsize,double mina,double maxa,double minp,double maxp,double(*f)(double),double x0,double x1){
    double xs[numpoints];
    double ys[numpoints];
    Descender::genPoints(numpoints,x0,x1,f,xs,ys);

    Fourier d(size,numpoints,xs,ys,mina,maxa,minp,maxp);

    cout << "Fourier Series Gradient Decent" <<
"\nSeries Size: " << size <<
"\nNum Points : " << numpoints <<
"\nIterations : " << iterations <<
"\nTolerance  : " << tolerance <<
"\nStep Size  : " << stepsize <<
"\nmin a      : " << mina <<
"\nmax a      : " << maxa <<
"\nmin p      : " << minp <<
"\nmax p      : " << maxp <<
"\nx0         : " << x0 <<
"\nx1         : " << x1 << "\n\n";

    cout << "Data x's:\n";
    for(int i = 0; i < numpoints; i++){
        cout << xs[i] << ", ";
    }
    cout << "\n\nData y's:\n";
    for(int i = 0; i < numpoints; i++){
        cout << ys[i] << ", ";
    }
    cout << "\n\nIterations:\n";

    double cost;
    for(int i = 0; i < iterations; i++){
        cost = d.g();
        if(i%10 == 0){
            cout << i << ": " << cost << "\n";
        }
        if(cost <= tolerance){
            cout << "\nConverged at epoch " << i << "\n";
            break;
        }
        d.iterate(stepsize);
    }
    cout << iterations << ": " << d.g() << "\n";

    cout << "\n\ncoefficients:\n";
    for(int n = 0; n < size; n++){
        cout << d.a[n] << ", ";
    }
    cout << "\n";

}

int main(int argc, char* argv[]){
    srand(time(NULL));

    const int TERMS = 100;
    const int ITERATIONS = 100;
    const double STEPSIZE = 0.1;
    const double TOLERANCE = 0.001;
    const unsigned int NUMPOINTS = 100;

    //Poly stuff
    const double PX0 = -1;
    const double PX1 = 0;
    const double MINK = -2;
    const double MAXK = 2;
    double (*POLYFUNCTION)(double) = polynomialfunct;
    PolynomialDescent(TERMS,NUMPOINTS,ITERATIONS,TOLERANCE,STEPSIZE,MINK,MAXK,POLYFUNCTION,PX0,PX1);

    //Fourier stuff
    const double FX0 = 0;
    const double FX1 = 10;
    const double MINA = -2;
    const double MAXA = 2;
    const double MINP = -2;
    const double MAXP = 2;
    double (*FOURIERFUNCTION)(double) = smoothfunct;
    FourierDescent(TERMS,NUMPOINTS,ITERATIONS,TOLERANCE,STEPSIZE,MINA,MAXA,MINP,MAXP,FOURIERFUNCTION,FX0,FX1);
}
