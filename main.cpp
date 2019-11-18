#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

struct Decender{
    unsigned int N;
    unsigned int D;
    double* P;
    double* a;
    double* d;
    double* y;

    Decender(unsigned int size, unsigned int numdata,double* datain, double* dataout, double mina, double maxa, double p0, double p1)
    :N(size),D(numdata),d(datain),y(dataout){
        a = new double[N];
        P = new double[N];
        for(int n = 0; n < N; n++){
            a[n] = ((double)rand()/RAND_MAX)*(maxa-mina) + mina;
            P[n] = p0 + ((1.0*n/(N-1)) * (p1 - p0));
        }
    }
    Decender(){}
    ~Decender(){}

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

    double dgdan(unsigned int n){
        double sum = 0;
        for(int i = 0; i < D; i++){
            sum += (2.0/D) * (k(d[i]) - y[i]) * sin(d[i] * P[n]);
        }
        return sum;
    }

    void iterate(double stepsize){
        double* newa = new double[N];
        for(int n = 0; n < N; n++){
            newa[n] = a[n] - stepsize*dgdan(n);
        }
        for(int n = 0; n < N; n++){
            a[n] = newa[n];
        }
        delete[] newa;
    }
};

void genPoints(uint length, double x0, double x1, double (*f)(double), double* d, double* y){
    for(int i = 0; i < length; i++){
        d[i] = x0 + ((1.0*i/(length-1)) * (x1 - x0));
        y[i] = (*f)(d[i]);
    }
}

void decend(unsigned int size,unsigned int numpoints,unsigned int iterations,double tolerance,double stepsize,double mina,double maxa,double minp,double maxp,double(*f)(double),double x0,double x1){
    double xs[numpoints];
    double ys[numpoints];
    genPoints(numpoints,x0,x1,f,xs,ys);

    Decender d(size,numpoints,xs,ys,mina,maxa,minp,maxp);

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

/*square wave with period = 2
f(x) = (floor(x)%2)*2 - 1

 1.  --  --
 0..1.2.3.4
-1.--  --
*/
double squarewave(double x){
    return 2*((int)x%2) - 1;
}

double smoothfunct(double x){
    return 0.5*sin(x) + 0.5*sin(x/2) + sin(-2/3.0*x);
}

int main(int argc, char* argv[]){
    srand(time(NULL));

    const int SIZE = 10;
    const int ITERATIONS = 1000;
    const double STEPSIZE = 0.05;
    const double TOLERANCE = 0.001;
    const unsigned int NUMPOINTS = 100;
    const double MINA = -2;
    const double MAXA = 2;
    const double MINP = -2;
    const double MAXP = 2;
    const double X0 = 0;
    const double X1 = 10;
    double (*FUNCTION)(double) = smoothfunct;

    decend(SIZE,NUMPOINTS,ITERATIONS,TOLERANCE,STEPSIZE,MINA,MAXA,MINP,MAXP,FUNCTION,X0,X1);

    return 0;
}
