#pragma once
#include <vector>
#include <ostream>
#include <cstdlib>
#include <cmath>
class band_matrix{
public:
    band_matrix(int n,int l);
    double get(int i,int j) const;
    void set(int i,int j,double v);
    void generate(int seed,double q);
    band_matrix T()const ;
    band_matrix operator*(const band_matrix &other)const;
    std::vector<double> operator*(const std::vector<double> &other)const;

protected:
    std::vector<std::vector<double>> bands;
    int n;
    int l;
};
struct linsys{
    band_matrix A;
    std::vector<double> b;
};
linsys gensys(int n,int l,double q,int seed);
std::vector<double> gen_vec(int n);
double operator*(const std::vector<double>&a,const std::vector<double>&b);
std::vector<double> operator-(const std::vector<double>&a,const std::vector<double>&b);
std::vector<double> operator+(const std::vector<double>&a,const std::vector<double>&b);
std::vector<double> operator*(double k,const std::vector<double>&b);
double norm(const band_matrix &A,const std::vector<double> &v);
double norm(const std::vector<double> &v);
std::ostream &operator<<(std::ostream &out,const std::vector<double> &vec);
