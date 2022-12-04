#include "generate.hpp"
#include <iostream>

int SOR(linsys & sys,std::vector<double>& x,double omega,int iters){
    for(int it = 0;it<iters;it++){
        std::vector<double> x0 = x;
        for(int i=0;i<sys.b.size();i++){
            double sum = 0.0;
            for(int j=0;j<sys.b.size();j++){
                if(i==j)continue;
                sum-=sys.A.get(i,j)*x[j];
            }
            sum+=sys.b[i];
            x[i] =x0[i]*(1-omega)+ sum*omega/sys.A.get(i,i);
        }
        if(norm(sys.A,x-x0)<1e-5)
            return it;
    }
    return iters;
}
const int N = 800;
const int L = 23;
int main(){
    auto sys = gensys(N,L,1.1,1234);
    auto sys2 = gensys(N,L,2,755);
    auto sys3 = gensys(N,L,10,234);
    for(int k=1;k<=7;k++){
        double omega = 2.0*k/8;
        std::vector<double> x(N,0);
        std::vector<double> x2(N,0);
        std::vector<double> x3(N,0);
        std::cout<<"omega = "<< omega<<" q = 1.1 n="<<SOR(sys,x, omega,100)<<std::endl;
        std::cout<<"omega = "<< omega<<" q = 2   n="<<SOR(sys2,x2, omega,100)<<std::endl;
        std::cout<<"omega = "<< omega<<" q = 10  n="<<SOR(sys3,x3,omega,100)<<std::endl;
    }
    
}
