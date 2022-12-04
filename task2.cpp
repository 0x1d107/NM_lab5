#include "generate.hpp"
#include <iostream>

int jacobi(linsys & sys,std::vector<double>& x,int iters){
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
const int N = 800;
const int L = 23;
int main(){
    auto sys = gensys(N,L,1.1,1234);
    auto sys2 = gensys(N,L,2,755);
    auto sys3 = gensys(N,L,10,234);
    
    std::vector<double> x(N,0);
    std::vector<double> x2(N,0);
    std::vector<double> x3(N,0);
    std::cout<<"q = 1.1 n="<<jacobi(sys,x,100)<<std::endl;
    std::cout<<"q = 2 n="<<jacobi(sys2,x2,100)<<std::endl;
    std::cout<<"q = 10 n="<<jacobi(sys3,x3,100)<<std::endl;
    
}
