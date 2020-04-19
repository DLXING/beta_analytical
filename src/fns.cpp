// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;

/// [[Rcpp::export]]
CharacterVector csample_char(const CharacterVector& x) {
    // same as sample in R, but for CharacterVector
    CharacterVector ret = 
        Rcpp::RcppArmadillo::sample(x, x.size(), false, NumericVector::create()) ;
    return ret ;
}
/// [[Rcpp::export]]
NumericMatrix table2(const CharacterVector& x, const CharacterVector& y){
    // same as table in R, but for CharacterVector
    CharacterVector xlevs =sort_unique(x);
    CharacterVector ylevs =sort_unique(y);
    IntegerVector xi = match(x, xlevs)-1;
    IntegerVector yi = match(y, ylevs)-1;
    NumericMatrix res(xlevs.size(),ylevs.size());
    for(int i=0; i<x.size(); i++)
        res(xi[i],yi[i]) +=1.0;
    rownames(res)=xlevs;
    colnames(res)=ylevs;
    return res;
}
// [[Rcpp::export]]
NumericMatrix commsimu(const DataFrame& ind){
    // function for the null model of Kraft et al. 2011. 
    //  i.e. shuffleing the species - subplot relationship.
    // @ind is a dataframe with each row a individual. 
    //  it has two columns sp and subplot. 
    
    CharacterVector swap = csample_char(ind("subplot"));
    NumericMatrix comp=table2(ind("sp"),swap);
    return comp;
}
// [[Rcpp::export]]
double Beta(const NumericMatrix& comp){
    // function for calculating beta. this function uses Chao & Chiu's (2016)
    //  formulas, but the results are identical to that from the 
    //  beta_obs R function.
    // @comp- S*M community matrix with either occurence or abundance data
    double M=comp.ncol();

    NumericVector p = rowSums(comp)/sum(comp);
    double Dg, Da;
    p=p[p>0];
    
    NumericVector pa = comp/sum(comp);
    pa=pa[pa>0];
    
    Dg=p.size();//zero-order gamma diversity, sensu Chao&Chiu 2016
    Da=pa.size()/M;//zero-order alpha diversity, sensu Chao&Chiu 2016
    double res= 1.0 - Da/Dg;//Beta_W=M/(M-1)*Beta_J

    return res;
}

// [[Rcpp::export]]
arma::mat do_simu(List ind_dat, int nsim=999){
    // function for calculating simulated beta. the function returns a matrix 
    //  with nrow equals to length of ind_dat, and ncol equals to nsim. 
    // @ind_dat- list of individual-level data. see mat_to_int for details.
    // @nsim- number of simulations.
    // 
    arma::mat bp(ind_dat.size(),nsim);
    
    for(int i=0; i<ind_dat.size(); i++){
        DataFrame ind_i(ind_dat[i]);
        for(int j=0; j<nsim; j++){
            NumericMatrix comp=commsimu(ind_i);
            bp(i,j) = Beta(comp);
        }
        if(i%20==1)
            Rcout<<i<<"\t"; //printf("%i\t",  i); //
    }
    printf("\n");
    return bp;
}
