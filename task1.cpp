#include <iostream>
#include "generate.hpp"
const int N =20;
const int L = 2;
int main(){
    linsys S1 = gensys(N,L,1.1,1234);
    linsys S2 = gensys(N,L,2,5678);
    linsys S3 = gensys(N,L,10,9876);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            std::cout<<S1.A.get(i,j)<< ' ';
        }
        std::cout<<std::endl;
    }
    std::cout <<" =\n"<<S1.b<<std::endl;
    
    
    
    return 0;
}
