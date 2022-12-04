#include "generate.hpp"
#include <iostream>


int jacobi(const linsys & sys,std::vector<double>& x,int iters){
    for(int it = 0;it<iters;it++){
        std::vector<double> x0 = x;
        for(int i=0;i<sys.b.size();i++){
            double sum = 0.0;
            for(int j=0;j<sys.b.size();j++){
                if(i==j)continue;
                sum-=sys.A.get(i,j)*x0[j];
            }
            sum+=sys.b[i];
            x[i] = sum/sys.A.get(i,i);
        }
        if(norm(sys.A,x-x0)<1e-5)
            return it;
    }
    return iters;
}
int CGM(linsys & sys,std::vector<double>& x,int iters){
    std::vector<double> p = sys.b - sys.A*x;
    for(int it=0;it<iters;it++){
        auto q = sys.A*p;
        auto r = sys.b - sys.A*x;
        auto alpha = - (r*r)/(p*q);
        x = x - alpha*p;
        auto r0 = r;
        r = r + alpha*q;
        if(norm(r)<1e-5)
            return it;
        double beta = (r*r)/(r0*r0);
        p = r+beta*p;
        

    }
    return iters;
}
int PCGM(linsys & sys,std::vector<double>& x,int m,int iters){
    std::vector<double> r0 = sys.b - sys.A*x;
    std::vector<double> R0(sys.b.size(),0); 
    jacobi({sys.A,r0},R0,m);
    std::vector<double> p = R0;

    for(int it=0;it<iters;it++){
        auto q = sys.A*p;

        auto alpha = - (R0*r0)/(p*q);
        x = x - alpha*p;
        auto r = r0 + alpha*q;
        if(norm(r)<1e-5)
            return it;
        std::vector<double> R(sys.b.size(),0); 
        jacobi({sys.A,r},R,m);
        double beta = (R*r)/(R0*r0);
        p = R+beta*p;
        R0 = R;
        r0 =r;

    }
    return iters;
}
const int N = 800;
const int L = 23;
int main(){
    auto sys = gensys(N,L,1.1,1234);
    auto sys2 = gensys(N,L,2,755);
    auto sys3 = gensys(N,L,10,234);
    
    std::vector<double> x(N,0);
    std::vector<double> x2(N,0);
    std::vector<double> x3(N,0);
    std::cout<<"q = 1.1 n="<<CGM(sys,x,100)<<std::endl;
    std::cout<<"q = 2 n="<<CGM(sys2,x2,100)<<std::endl;
    std::cout<<"q = 10 n="<<CGM(sys3,x3,100)<<std::endl;
    for(int m=1;m<10;m++){
        
        std::vector<double> xp(N,0);
        std::cout<<"m = "<<m<<" q = 10 n="<<PCGM(sys3,xp,m,100)<<std::endl;
    }
    
}
