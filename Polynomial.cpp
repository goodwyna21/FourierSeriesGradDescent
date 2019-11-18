#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

struct Decender{
    unsigned int N;
    unsigned int D;
    double* k;
    double* d;
    double* y;

    Decender(unsigned int degree, unsigned int numdata,double* datain, double* dataout, double mink, double maxk)
    :N(degree),D(numdata),d(datain),y(dataout){
        k = new double[N+1];
        for(int n = 0; n <= N; n++){
            k[n] = ((double)rand()/RAND_MAX)*(maxk-mink) + mink;
        }
    }
    Decender(){}
    ~Decender(){}

    double h(double x){
        double sum = 0;
        for(int n = 0; n <= N; n++){
            sum += k[n] * pow(x,n);
        }
        return sum;
    }

    double g(){
        double sum = 0;
        for(int i = 0; i < D; i++){
            sum += (1.0/D)*pow(h(d[i]) - y[i], 2);
        }
        return sum;
    }

    double dgdkn(unsigned int n){
        double sum = 0;
        for(int i = 0; i < D; i++){
            sum += (2.0/D) * (h(d[i]) - y[i]) * pow(d[i],n);
        }
        return sum;
    }

    void iterate(double stepsize){
        double* newk = new double[N];
        for(int n = 0; n <= N; n++){
            newk[n] = k[n] - stepsize*dgdkn(n);
        }
        for(int n = 0; n <= N; n++){
            k[n] = newk[n];
        }
    }
};

void genPoints(uint length, double x0, double x1, double (*f)(double), double* d, double* y){
    for(int i = 0; i < length; i++){
        d[i] = x0 + ((1.0*i/(length-1)) * (x1 - x0));
        y[i] = (*f)(d[i]);
    }
}

void decend(unsigned int size,unsigned int numpoints,unsigned int iterations,double tolerance,double stepsize,double mina,double maxa,double(*f)(double),double x0,double x1){
    double xs[numpoints];
    double ys[numpoints];
    genPoints(numpoints,x0,x1,f,xs,ys);

    Decender d(size,numpoints,xs,ys,mina,maxa);

    cout << "Fourier Series Gradient Decent" <<
"\nSeries Size: " << size <<
"\nNum Points : " << numpoints <<
"\nIterations : " << iterations <<
"\nTolerance  : " << tolerance <<
"\nStep Size  : " << stepsize <<
"\nmin a      : " << mina <<
"\nmax a      : " << maxa <<
"\nx0         : " << x0 <<
"\nx1         : " << x1 << "\n\n";

    cout << "Initial coefficients:\n";
    for(int n = 0; n < size+1; n++){
        cout << d.k[n] << ", ";
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
        cout << d.k[n] << ", ";
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

double polynomial(double x){
    return -0.8*pow(x,3) - 2*pow(x,2) - x + 0.6;
}

int main(int argc, char* argv[]){
    srand(time(NULL));

    const int DEGREE = 100;
    const int ITERATIONS = 1000;
    const double STEPSIZE = 0.1;
    const double TOLERANCE = 0.001;
    const unsigned int NUMPOINTS = 100;
    const double X0 = -1;
    const double X1 = 0;
    const double MINK = -2;
    const double MAXK = 2;
    double (*FUNCTION)(double) = polynomial;

    decend(DEGREE,NUMPOINTS,ITERATIONS,TOLERANCE,STEPSIZE,MINK,MAXK,FUNCTION,X0,X1);

    return 0;
}
