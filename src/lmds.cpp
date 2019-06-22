#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <list>

using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

float cc = 0;
float t = 0;
const float epsilon = 0.00001;
float stepsize = 0.1;
int n = 0; //number of rows in X
float Mka = 0;
float s0 = 0;
float s1 = 0;

//instantiate all the necessary matrices
mat Do;
mat Daux;
mat Inb;
mat Inb1;
mat D1mulam;
mat D1mu;
mat k_v;
mat Dnu;
mat Dnulam;
mat D1mu2;
mat D1mulam2;
mat M;
mat D1_Inb;
mat D1_rk;
mat Nk;

mat Grad;
mat X0;
mat D1;
mat pc1;
mat pc2;
mat pca;

mat getDistanceMatrix(mat X){
    mat dm(n,n,fill::zeros); vec X2;
    X2 = vectorise(sum(pow(X,2.0),1));
    dm.each_row() += trans(X2);
    dm.each_col() += X2;
    dm-= 2 * X * trans(X);
    dm.elem(find(dm<0.0)).zeros();
    return sqrt(dm);
};

mat comparisonMatrix(mat X1, mat X2){
    mat boolMatrix(n,n,fill::zeros);
    for(int i = 0; i<n;i++){
        for(int j = 0; j<n; j++){
            boolMatrix(i,j) = X1(i,j) - X2(0,i) > 1e-10 ? 0 : 1;
        }
    }
    return boolMatrix;
};

mat comparisonMatrix2(umat X1, mat X2){
    mat boolMatrix(n,n,fill::zeros);
    for(int i = 0; i<n;i++){
        for(int j = 0; j<n; j++){
            boolMatrix(i,j) = (float) X1(i,j) - X2(0,i) > 1e-10 ? 0 : 1;
        }
    }
    return boolMatrix;
};

mat parallelMax(mat X){
    //This function is made specific with respect to the fact that
    //the input matrix consists of only 1's and 0's ideally.
    mat result(n,n,fill::zeros);
    result = X + trans(X);
    result.elem(find(result>1.1)) -= 1.0;
    return result;
};

mat specificExponent(mat X1,mat X2, float exponentValue){
    mat result(n,n);
    if(exponentValue!=0.0){
        result = pow(X1 % X2,exponentValue);
        result.elem(find(X1<0.9)) *=0;
        result.diag() *= 0;
    }
    else{
        result = X1;
        result.diag() *= 0.0;
    }
    return result;
};


mat rankRows(arma::mat X){
//different rank method that returns first appearance.
    int n = X.n_rows;
    mat rankedMatrix(n,n,fill::zeros);
    vec temp;
    for(int i = 0; i<n; i++){
        temp = vectorise(sort(X.row(i),"ascend"));
        for(int j = 0; j<n; j++){
            for(int k = 0; k<n; k++){
                if(temp(k)==X(i,j)){
                    rankedMatrix(i,j) = k;
                    break;
                }
            }
        }
    }
    return rankedMatrix;
}

mat rankRows2(arma::mat X){
//This does not do anything for tied values. It just continues
//to order with no duplicates or averages.
    umat rankedMatrix(n,n,fill::zeros);
    uvec temp;
    for(int i = 0; i<n; i++){
        temp = vectorise(sort_index(X.row(i),"ascend"));
        rankedMatrix.row(i) += temp.t();
    }
    mat rankedMat = conv_to<mat>::from(rankedMatrix);
    return rankedMat;
};

float twoNorm(mat X){
    return sqrt(accu(pow(X,2)));
};

mat mds_c(mat X,int dim = 2){
    int n = X.n_rows;
    int m = X.n_cols;
    mat X2 = X;
    X2.each_row() -= mean(X,0);
    X2.each_row() /= stddev(X,0);
    mat D(n,n,fill::zeros);
    D = getDistanceMatrix(X2);
    D = square(D);
    mat C(n,n,fill::zeros);
    C.diag().ones();
    C -= 1.0/float(n);
    mat B(n,n,fill::zeros);
    B = -0.5 * C * D * C;
    vec eigval;mat eigvec;
    eig_sym(eigval,eigvec,B);
    mat G2(n,n,fill::zeros);
    G2.diag() += sqrt(eigval);
    mat Xk = eigvec(span(0,n-1),span(m-dim+1,m));
    Xk = Xk * G2(span(m-dim+1,m),span(m-dim+1,m));
    return Xk;
};




// [[Rcpp::export]]
Rcpp::List localMDS_c(arma::mat X, Rcpp::Nullable<arma::mat> X1_temp = R_NilValue, int random_start = 0, float k = 6.0, float d = 2.0, float lambda = 1.0,
              float mu = 1.0, float nu = 0.0, float tau = 1.0, int niter = 500){
    std::list<float> stress;
    Rcpp::List ret;
    Mka = 0;
    n = X.n_rows;
    Do = getDistanceMatrix(X);
    mat E(n,d,fill::ones);
    Daux = sort(Do,"ascend",0);
    Daux = Daux.row(k);
    Inb = comparisonMatrix(Do,Daux);
    k_v = sum(Inb,1);
    k = (accu(k_v)-n)/n;
    mat Inbsum(n,n);
    Inbsum.each_col() = vectorise(k_v);
    Inb1 = parallelMax(Inb);
    Dnu = specificExponent(Inb1,Do,nu);
    Dnulam = specificExponent(Inb1,Do,(nu+1/lambda));
    Dnu.diag().zeros();
    Dnulam.diag().zeros();
    cc = (accu(Inb1)-n)/n/n*median(Dnulam.elem( find(Dnulam!=0)));
    t = tau * cc;
    mat X1;
    if((X1_temp.isNull())&&(random_start==0)){
        X1 = randn(n,d);
    }
    else if((X1_temp.isNull())&&(random_start==1)){
        X1 = mds_c(X,d);
    }
    else if(X1_temp.isNotNull()){
        NumericMatrix X2_temp(X1_temp.get());
        X1 = as<mat>(X2_temp);
    }
    D1 = getDistanceMatrix(X1);
    X1 = X1 * twoNorm(Do)/twoNorm(D1);
    s1 = std::numeric_limits<float>::infinity();
    s0 = 0;
    stepsize = 0.1;
    int i = 0;
    while((stepsize> epsilon) && (i < niter)){
        if((s1>=s0)&&(i>1)){
            stepsize = stepsize * 0.5;
            X1 = X0 - stepsize * Grad;
        }
        else{
             stepsize = stepsize * 1.05;
             X0 = X1;
             D1mu2 = pow(D1,mu-2.0);
             D1mu2.diag().zeros();
             D1mulam2 = pow(D1,(mu+1/lambda - 2.0));
             D1mulam2.diag().zeros();
             M =  Dnu % D1mulam2 - D1mu2 % (Dnulam + t*-1*(Inb1-1));
             Grad = X0%(M*E) - M*X0;
             Grad = (twoNorm(X0)/twoNorm(Grad))*Grad;
             X1 = X0 - stepsize*Grad;
        }
        i = i+1;
        s0 = s1;
        D1 = getDistanceMatrix(X1);
        D1mulam = pow(D1,(mu+1/lambda));
        D1mulam.diag().zeros();
        D1mu = pow(D1,mu);
        D1mu.diag().zeros();
        //only check other scenarios if one is not witnessed...
        if((mu+1/lambda)==0)
        {
            D1.diag().ones();
            s1 = accu(Dnu%log(D1))-accu((D1mu-1)%Dnulam)/mu - t*accu((D1mu-1)%(1-Inb1))/mu;
        }
        else{
            if(mu==0)
            {
                D1.diag().ones();
                s1 = accu(Dnu%(D1mulam-1))/(mu+1/lambda) - accu(log(D1)%Dnulam)-t*accu(log(D1)%(1-Inb1));
            }
            else{
                if((mu!=0)&&((mu+1/lambda)!=0))
                {
                    s1 = accu(Dnu%(D1mulam-1))/(mu+1/lambda) - accu((D1mu-1)%Dnulam)/mu - t*accu((D1mu-1)%(1-Inb1))/mu;
                }
            }
        }
        if((i+1)%100==0){
            stress.push_back(s1);
            D1_rk = rankRows2(D1);
            D1_Inb = comparisonMatrix(D1_rk,Inbsum);
            Nk = (sum(D1_Inb%Inb,1)-1)/(k_v-1);
            Mka = accu(Nk)/n - k/n;
            std::cout<<"niter="<<(i+1)<<" Stress="<<s1<<" Mka_adj="<<Mka<<endl;
        }
    }
    //Want to return the first two principal components
    pca = princomp(X1);
    pc1 = pca.col(0);
    pc2 = pca.col(1);
    ret["X"] = X1;
    ret["pc1"] = pc1;
    ret["pc2"] = pc2;
    ret["stress"] = stress;
    ret["Mka"] = Mka;
    return ret;
}
