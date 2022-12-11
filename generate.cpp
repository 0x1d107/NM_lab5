#include "generate.hpp"
#include <cmath>
band_matrix::band_matrix(int n,int l){
    this->n = n;
    this->l = l;
    bands = std::vector<std::vector<double>>(n,std::vector<double>(2*l+1,0.0));
    
}
double band_matrix::get(int i,int j)const{
    int x = j - i + l;
    int y = i;
    if(y>=0&&y<n&&x>=0&&x<2*l+1)
        return bands[y][x];
    return 0;
}
void band_matrix::set(int i,int j, double v){
    int x = j - i + l;
    int y = i;
    if(y>=0&&y<n&&x>=0&&x<2*l+1)
        bands[y][x]=v;
    
}
double rnd(){
    return rand()*2.0/RAND_MAX - 1.0;
}
void band_matrix::generate(int seed,double q){
    srand(seed);
    for(int i=0;i<n;i++){
        double sum =0.0;
        for(int j =0;j<2*l+1;j++){
            if(j==l)
                continue;
            double v = rnd();
            bands[i][j] = v;
            sum+=fabs(v);
        }
        bands[i][l] = q*sum;
    }

    
}
band_matrix band_matrix::operator*(const band_matrix &other)const {
    band_matrix result(n,2*l);
    for(int i=0;i<n;i++){
        for(int j=0;j<other.n;j++){
            double sum=0.0;
            for(int k=0;k<n;k++)
                sum+=get(i,k)*other.get(k,j);
            result.set(i,j,sum);
        }
    }
    return result;
    
}
std::vector<double> band_matrix::operator*(const std::vector<double> &other)const {
    std::vector<double> result(other.size());
    for(int i=0;i<n;i++){
        double sum=0.0;
        for(int k=0;k<n;k++)
            sum+=get(i,k)*other[k];
        result[i] = sum;
        
    }
    return result;
}
band_matrix band_matrix::T() const{
    band_matrix result(n,l);
    for(int i =0;i<n;i++){
        for(int j=0;j<n;j++){
            result.set(i,j,get(j,i));
        }
    }
    return result;
}
std::vector<double> gen_vec(int n){
    std::vector<double> res;
    for(int i=0;i<n;i++){
        res.push_back(rnd());
    }
    return res;
}
std::ostream &operator<<(std::ostream &out,const std::vector<double> &vec){
    out<<'[';
    for(int i=0;i<vec.size();i++){
        out<<vec[i];
        if(i !=vec.size()-1) 
            out<<',';
    }
    out<<']';
    return out;
}
linsys gensys(int n,int l,double q,int seed){
    linsys s = {{n,l},std::vector<double>()};
    s.A.generate(seed,q);
    auto x = gen_vec(n);
    s.b = s.A.T()*(s.A * x);
    s.A = s.A.T()*s.A;
    for(int i =0;i<n;i++){
        double sum = 0;
        for(int j=(i>=l?i-l:0);j<n && j<=i+l;j++){
            if(i==j)
                continue;
            sum+=fabs(s.A.get(i,j));
        }
        if(fabs(s.A.get(i,i))< sum)
            s.A.set(i,i,sum);

    }
    return s;
}

double operator*(const std::vector<double>&a,const std::vector<double>&b){
    double sum = 0;
    for(int i =0;i<a.size();i++){
        sum+=a[i]*b[i];
    }
    return sum;
}
std::vector<double> operator-(const std::vector<double>&a,const std::vector<double>&b){
    std::vector<double> c(a.size());
    for(int i=0;i<a.size();i++){
        c[i] = a[i] - b[i];
    }
    return c;
}
std::vector<double> operator+(const std::vector<double>&a,const std::vector<double>&b){
    std::vector<double> c(a.size());
    for(int i=0;i<a.size();i++){
        c[i] = a[i] + b[i];
    }
    return c;
}
std::vector<double> operator*(double k,const std::vector<double>&a){
    std::vector<double> c(a.size());
    for(int i=0;i<a.size();i++){
        c[i] = a[i] *k;
    }
    return c;
}
double norm(const band_matrix &A,const std::vector<double> &v){
    return sqrt((A*v)*v);
}
double norm(const std::vector<double> &v){
    return sqrt(v*v);
}
