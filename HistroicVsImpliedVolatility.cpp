#include <iostream>
#include "mvector.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>

// payoff function declaration
double Payoff_European_Call(double y, double X);

// payoff function declaration
double Payoff_European_Put(double y, double X);

// A function declaration
double A(double x, double r, double sigma, double T, double k);

// B function declaration
double B(double x, double y, double sigma, double T, double k);

// base integral class
class Integral
{
    public:
        virtual double operator()(const double& y)=0;
        virtual double residual(const double& x, double exact)=0;
};

// integrate function declaration
double integrate(Integral &f, double a, double b, int n);

//derived test integral class
class Test_Integral: public Integral
{
    public:
        Test_Integral(){X=10; r=0.05; sigma=0.1;T=0.5;}
        virtual double operator()(const double& y) 
        {
            return X*exp(T*(y-r*sigma));
        }
        void SetX(double X_){X=X_;}
        void Setr(double r_){r=r_;}
        void Setsigma(double sigma_){sigma=sigma_;}
        void SetT(double T_){T=T_;}
        private:
        double X, r, sigma, T;
};

// derived european option class
class European_Call: public Integral
{
    public:
        European_Call(){X=10; r=0.05; sigma=0.1;T=0.5; y=0.5; S_0=10;}
        virtual double operator()(const double& y) 
        {
            double k=(2*r)/pow(sigma,2)-1, x=std::log(S_0/X);
            return A(x,r,sigma,T,k)*B(x,y,sigma,T,k)*Payoff_European_Call(y,X);
        }
        virtual double residual(const double& x, double exact) 
        {
            sigma=x;
            double a = -10.0;
            double b = 10.0;
            int n = 1000; // must be even
            return integrate(*this, a, b, n)-exact;
        }
        void SetX(double X_){X=X_;}
        void Setr(double r_){r=r_;}
        void Setsigma(double sigma_){sigma=sigma_;}
        void SetT(double T_){T=T_;}
        void Sety(double y_){y=y_;}
        void SetS_0(double S_0_){S_0=S_0_;}
        private:
        double X, r, sigma, T, y, S_0;
};

// derived european option class
class European_Put: public Integral
{
    public:
        European_Put(){X=10; r=0.05; sigma=0.1;T=0.5; y=0.5; S_0=10;}
        virtual double operator()(const double& y) 
        {
            double k=(2*r)/pow(sigma,2)-1, x=std::log(S_0/X);
            return A(x,r,sigma,T,k)*B(x,y,sigma,T,k)*Payoff_European_Put(y,X);
        }
        virtual double residual(const double& x, double exact) 
        {
            sigma=x;
            double a = -10.0;
            double b = 10.0;
            int n = 1000; // must be even
            return integrate(*this, a, b, n)-exact;
        }
        void SetX(double X_){X=X_;}
        void Setr(double r_){r=r_;}
        void Setsigma(double sigma_){sigma=sigma_;}
        void SetT(double T_){T=T_;}
        void Sety(double y_){y=y_;}
        void SetS_0(double S_0_){S_0=S_0_;}
        private:
        double X, r, sigma, T, y, S_0;
};

// derived quadratic test function class
class Quadratic_f: public Integral
{
    public:
        Quadratic_f(){}
        virtual double residual(const double& x, double exact)
        {
            exact=0;
            return x*x-2.0-exact;
        } 
        virtual double operator()(const double& x) 
        {
            return x*x-2.0;
        }
};

// estimates yearly variance of a stock following a log normal random walk
double Variance_Estimator(MVector v){
    MVector w;
    for (int i=1;i<v.size();i++)
    {
        w.push_back(std::log(v[i])-std::log(v[i-1]));
    }
    double sum=0, sum_of_squares=0;
    for (int i=0;i<w.size();i++)
    {
        sum+=w[i];
        sum_of_squares+=w[i]*w[i];
    }
    double var;
    var=(sum_of_squares/(w.size()-1))-((sum*sum)/(w.size()*w.size()-w.size()));
    return std::sqrt(252*var);
}

// function to implement Simpsons rule from y=a to y=b with n steps
double integrate(Integral &f, double a, double b, int n)
{
    if (n%2!=0 || a>b) {std::cout<<"Error"<<std::endl; return 1;}
    double sum=0, h=(b-a)/n;
    sum+=f(a)+f(b);
    for (int j=1; j<=n/2-1;j++)
    {
        sum+=2*f(a+2*j*h)+4*f(a+(2*j-1)*h);
    }
    sum=(h/3)*(sum+4*f(a+(n-1)*h));
    return sum;
}


double Payoff_European_Call(double y, double X)
{
    return std::max(X*(exp(y)-1),0.0);
}

double Payoff_European_Put(double y, double X)
{
    return std::max(X*(1-exp(y)),0.0);
}

double A(double x, double r, double sigma, double T, double k)
{
    return exp(-0.5*k*x-0.125*pow(sigma,2)*pow(k,2)*T-r*T)/std::sqrt(2*pow(sigma,2)*M_PI*T);
}

double B(double x, double y, double sigma, double T, double k)
{
    return exp(-pow(x-y,2)/(2*pow(sigma,2)*T)+0.5*k*y);
}

// secant method tester function

double newton_secant_f(double x) // returns f(x)
{
return x*x-2.0;
}

// secant iterator
bool NewtonSecant(Integral &f,double &result, double x0, double x1, double tol, int
maxiter, double exact)
{
    int i;
    double x;
    for (i=0; i<maxiter; i++)
    {
        x = x1 - f.residual(x1,exact)*(x1-x0)/
        (f.residual(x1,exact)-f.residual(x0,exact));
        if (std::abs(f.residual(x,exact))<tol) break;
        x0=x1;
        x1=x;
    }
    result = x;
    if (i==maxiter)
        return false;
    else return true;
}


int main()
{
    /*
    {
        std::ifstream input("ezj_open_prices.txt");
        MVector column1;
        if (!input)
        {
            std::cout << "Could not open file for reading" << std::endl;
            return 1;
        }
        // for (int i=0; i<5; i++)
        while (input)
        {
            double x1=0;
            if (input >> x1)
            {
                column1.push_back(x1);
            }
        }
        input.close();
    std::cout<<column1.size()<<std::endl;
    std::cout<<Variance_Estimator(column1)<<std::endl;
    }
    */
   
    /* 2.3.1 test a single case:
    {
    std::cout.precision(16);
    double X=10, r=0.05, sigma=0.1, T=0.5;
    Test_Integral f;
    int n=2;
    std::cout<<integrate(f,0,1,n)<<std::endl;
    std::cout<<0.5*(f(0)+4*f(0.5)+f(1))/3<<std::endl;
    double I=((X*exp(-r*sigma*T))*(exp(T)-1))/T;
    std::cout<<I-integrate(f,0,1,n)<<std::endl;
    }
    */

    /* 2.3.1 error analysis:
   {
    double X=10, r=0.05, sigma=0.1, T=0.5;
    double I=((X*exp(-r*sigma*T))*(exp(T)-1))/T;
    Test_Integral f;
    int n=2;
    for (int i=0;i<20;i++)
    {
        double h=1.0/n;
        double error_bound=(pow(T,4)*pow(h,4)*X*exp(T*(1-r*sigma)))/180;
        std::cout.precision(16);
        // std::cout<<integrate(f,0,1,n)<<std::endl;
        // std::cout<<I<<std::endl;
        // std::cout<<std::abs(I-integrate(f,0,1,n))<<std::endl;
        std::cout<<error_bound<<std::endl;
        n*=2;
    }
    }
    */

    /* multiple T
    {
        double X=10, r=0.05, sigma=0.1;
        int n=10e+6;
        Test_Integral f;
        std::cout.precision(16);
        for (int i=0; i<5; i++)
        {
            f.SetT(i/4.0);
            std::cout<<integrate(f, 0,1,n)<<std::endl;
            std::cout<<((X*exp(-r*sigma*(i/4.0)))*(exp(i/4.0)-1))/(i/4.0)<<std::endl;
        }

    }
    */
   /* Testing A,B and payoff
   {
    double X=10, r=0.05, sigma=0.1, T=0.5, y=0.5, S_0=10;
    double k=(2*r)/pow(sigma,2)-1, x=std::log(S_0/X);
    std::cout<<payoff(y,X)<<std::endl;
    std::cout<<A(x,r,sigma,T,k)<<std::endl ;
    std::cout<<B(x,y,sigma,T,k)<<std::endl ;
   }
   */

   /* Testing European_Call for different n
   {
    European_Call f;
    European_Put g;
    g.Setr(0.06);
    g.SetS_0(97);
    g.Setsigma(0.2);
    g.SetT(0.75);
    g.SetX(100);
    f.Setr(0.06);
    f.SetS_0(97);
    f.Setsigma(0.2);
    f.SetT(0.75);
    f.SetX(100);
    double exact_value=7.369373e+0;
    for(int n=200;n<=200;n+=200)
    std::cout<<integrate(f,-10.0,10.0,n)<<" "<<integrate(g,-10,10,n)<<std::endl;
   }
   */

   /* newton secant quadratic tester
   {
    Quadratic_f f;
    double guess=-1.5;
    std::cout<<NewtonSecant(f,guess,guess,guess+0.5,1e-6,2,0)<<std::endl;
    std::cout<<guess<<std::endl;
   }
   */ 

   /* newton secant call tester
   {
    European_Call f;
    double guess=0.5;
    double exact=1.425;
    f.Setr(0.002);
    f.SetS_0(30.23);
    f.SetT(0.1205);
    f.SetX(29);
    std::cout<<NewtonSecant(f, guess, guess, guess+0.1, 1e-8, 5, exact)<<std::endl;
    std::cout<<guess<<std::endl;
    
    for(int i=1;i<=30;i++)
    {
        f.Setsigma(i/100.0);
        std::cout<<integrate(f,-10, 10, 10000)-exact<<std::endl;
    }
   }
    */

   /* newton secant call tester with X variation
   {
    European_Call f;
    double guess=0.5;
    double exact=1.425;
    f.Setr(0.002);
    f.SetS_0(30.23);
    f.SetT(0.1205);
    f.SetX(29);
    std::cout<<NewtonSecant(f, guess, guess, guess+0.1, 1e-8, 5, exact)<<std::endl;
    std::cout<<guess<<std::endl;
    for(int i=0;i<=100;i++)
    {
        f.SetX(29+i/10.0);
        NewtonSecant(f, guess, guess, guess+0.1, 1e-8, 10, exact);
        std::cout<<guess<<std::endl;
    }
   }
   */

   /* newton secant call finder 
   {
    European_Call f;
    double guess=0.5;
    double exact=0.05;
    f.Setr(0.05);
    f.SetS_0(1.0);
    f.SetT(0.5);
    f.SetX(1);
    std::cout<<NewtonSecant(f, guess, guess, guess+0.1, 1e-8, 5, exact)<<std::endl;
    std::cout<<guess<<std::endl;
   }
    */



   /* amazon implied
   {
    European_Call f;
    double guess=2;
    double exact=(104.25+105.8)/2.0;
    f.Setr(0.04);
    f.SetS_0(104.55);
    f.SetT(1.0/365.0);
    f.SetX(130);
    std::cout<<NewtonSecant(f, guess, guess, guess+0.1, 1e-8, 10, exact)<<std::endl;
    std::cout<<guess<<std::endl;
   }
    */
   // error investigation
   {
    European_Call f;
    double guess=0.5;
    double exact=1.425;
    f.Setr(0.002);
    f.SetS_0(30.23);
    f.SetT(0.1205);
    f.SetX(29);
    std::cout<<NewtonSecant(f, guess, guess, guess+0.1, 1e-8, 5, exact)<<std::endl;
    std::cout<<guess<<std::endl;
    double limit=706;
    std::cout<<integrate(f,-limit, limit, limit*100)-exact<<std::endl;
   }

    return 0;
}