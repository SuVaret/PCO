#include <Rcpp.h>
#include <algorithm>
#include <RcppEigen.h>
#include <map>


using namespace Rcpp;
//[[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::Map;
using Eigen::VectorXd;




class criterion_1D{
public:
  Eigen::VectorXd xi;
  
  double hmin;
protected :
  int n;
  double n2;
  double hmin2;
  
public :
  criterion_1D(Eigen::VectorXd xi = Eigen::VectorXd::Zero(1)) {
    this->xi = xi;
    this->n = xi.size();
    this->n2 = n * n;
  }
  
public:
  virtual Eigen::VectorXd compute(Eigen::ArrayXd H) = 0;
};

class exact_crit_1D : public criterion_1D {
protected:
  Eigen::ArrayXd u; // u is supposed to be already squared
  
  
public :
  exact_crit_1D(Eigen::VectorXd xi):criterion_1D(xi){
    int size = n * (n - 1) / 2;
    this->u = Eigen::ArrayXd::Zero(size);
    
  }
  
protected :
  
  
  // Delta = 0 càd i!=j
  void outer_diff_square_1D(){
    // Computation of (oi - oj)^2 for all i and j s.t. j > i
    int pos_u = 0;
    int n_u = n - 1;
    for (int i = 0; i < (n - 1); i++){
      Rcpp::checkUserInterrupt();
      u.segment(pos_u, n_u) = (xi.segment(i + 1, n_u).array() - xi(i)).square();
      pos_u += n_u;
      n_u -= 1;
    }
  }
  
  
  
  void outer_diff_abs_1D(){
    // Computation of ABS(oi - oj) for all i and j s.t. j > i
    int pos_u = 0;
    int n_u = n - 1;
    for (int i = 0; i < (n - 1); i++){
      Rcpp::checkUserInterrupt();
      u.segment(pos_u, n_u) = (xi.segment(i + 1, n_u).array() - xi(i)).abs();
      pos_u += n_u;
      n_u -= 1;
    }
    std::sort(u.data(), u.data() + u.size());
  }
  
};




class GK_exact_crit_1D : public exact_crit_1D{
private :
  Eigen::VectorXd k;
  
public :
  GK_exact_crit_1D(Eigen::VectorXd xi) : exact_crit_1D(xi){
    this->hmin = M_1_SQRT_2PI / double(n);
    this->hmin2 = hmin * hmin;
    outer_diff_square_1D();
    this->k = Eigen::VectorXd::Zero(u.size());
  }
  
  GK_exact_crit_1D(Eigen::VectorXd xi, double hmin) : exact_crit_1D(xi){
    this->hmin = hmin;
    this->hmin2 = hmin * hmin;
    outer_diff_square_1D();
    this->k = Eigen::VectorXd::Zero(u.size());
  }
  
public:
  Eigen::VectorXd compute(Eigen::ArrayXd H){
    Eigen::ArrayXd H2_p_hmin2 = H.square() + hmin2;
    Eigen::VectorXd pen = M_SQRT_2dPI * (double(n) * H2_p_hmin2.sqrt()).inverse();
    
    Eigen::VectorXd loss(H.size());
    Eigen::VectorXd c1 = -(M_LN2 + M_LN_SQRT_PI + H.log());
    Eigen::VectorXd c2 = -(M_LN_SQRT_2PI + 0.5 * H2_p_hmin2.log());
    
    for (int no_h = 0; no_h < H.size(); no_h++){
      Rcpp::checkUserInterrupt();
      double h = H(no_h);
      double h2 = h * h; // compute the square is faster than H2 access
      
      double ch_1 = -0.25 / h2;
      k = (c1(no_h) + ch_1 * u).exp();
      double sk_1 = 2.0 * k.sum() + n / (2.0 * h * M_SQRT_PI);
      
      
      double h2_p_hm2 = h2 + hmin2; // sum is faster than H2_p_hmin2 access
      double ch_2 = -0.5 / h2_p_hm2;
      k = (c2(no_h) + ch_2 * u).exp();
      double sk_2 = 2.0 * k.sum() + n / (M_SQRT2 * M_SQRT_PI * std::sqrt(h2_p_hm2));
      
      double somme =  sk_1 - 2.0 * sk_2;
      loss(no_h) = somme;
    }
    loss = loss / n2;
    Eigen::VectorXd crit = loss + pen;
    return crit;
  }
  
  
};




class EK_exact_crit_1D : public exact_crit_1D{
public :
  EK_exact_crit_1D(Eigen::VectorXd xi) : exact_crit_1D(xi){
    this->hmin = 0.75 / double(n);
    this->hmin2 = hmin * hmin;
    outer_diff_abs_1D();
  }

  EK_exact_crit_1D(Eigen::VectorXd xi, double hmin) : exact_crit_1D(xi){
    this->hmin = hmin;
    this->hmin2 = hmin * hmin;
    outer_diff_abs_1D();
  }

public:
  Eigen::VectorXd compute(Eigen::ArrayXd H){
    Eigen::VectorXd pen = 1.5 * (H.square() - hmin2 / 5.0) * H.pow(-3) / double(n);

    Eigen::VectorXd loss(H.size());
    for (int no_h = 0; no_h < H.size(); no_h++){
      Rcpp::checkUserInterrupt();
      double h = H(no_h);
      // index for loop on u = (oi - oj)
      int u_ind = 0;
      double itu = u(u_ind);
      // multiplicative constant
      double c_0 = 0.0375 / h; // 0.0375 = 2 * 3 / 160
      double r = hmin / h;
      double r2 = r * r;
      double r3 = r2 * r;
      while ((itu <  h - hmin) & (u_ind < u.size())){// itu is in [0; h - hmin[
        // (Kh * Kh) (itu)
        double z = itu / h;
        double z2 = z * z;
        double z3 = z * z2;
        double z5 = z3 * z2;
        double k1 = c_0 * (32 - 40 * z2 + 20 * z3 - z5);

        // (Kh * Khmin) (itu)
        double k2 = 0.3 * (5 - r2 - 5 * z2) / h; // 0.3 = 2*3/20

        loss(no_h) += k1 - 2 * k2;
        u_ind++;
        itu = u(u_ind);
      }
      while ((itu <  h + hmin) & (u_ind < u.size())){// itu is in [h - hmin; h + hmin[
        // (Kh * Kh) (itu)
        double z = itu / h;
        double z2 = z * z;
        double z3 = z * z2;
        double z5 = z3 * z2;
        double k1 = c_0 * (32 - 40 * z2 + 20 * z3 - z5);

        // (Kh * Khmin) (itu)
        double p3 = std::pow(1 + r - z, 3);
        double p1 = 4 * (1 + r2) - 3 * (r * (z + 4) + z) - z2;
        double k2 = - c_0 * p3 * p1 / r3;

        loss(no_h) += k1 - 2 * k2;
        u_ind++;
        itu = u(u_ind);
      }
      while ((itu <  2 * h) & (u_ind < u.size())){// itu is in [h - hmin; h + hmin[
        // (Kh * Kh) (itu)
        double z = itu / h;
        double z2 = z * z;
        double z3 = z * z2;
        double z5 = z3 * z2;
        double k1 = c_0 * (32 - 40 * z2 + 20 * z3 - z5);

        // (Kh * Khmin) (itu) = 0

        loss(no_h) += k1;
        u_ind++;
        itu = u(u_ind);
      }
      // adds the terms s.t. i = j
      double k1_0 = n * c_0 * 16;
      double k2_0 = n * 0.15 * (5 - r2) / h; // 0.15 = 3 / 20
      loss(no_h) += k1_0 - 2 * k2_0;
      loss(no_h) /= n2;
    }
    Eigen::VectorXd crit = loss + pen;
    return(crit);
  }


};

class BK_exact_crit_1D : public exact_crit_1D{
public :
  BK_exact_crit_1D(Eigen::VectorXd xi) : exact_crit_1D(xi){
    this->hmin = 15.0 / (16.0 * double(n));
    this->hmin2 = hmin * hmin;
    outer_diff_abs_1D();
  }

  BK_exact_crit_1D(Eigen::VectorXd xi, double hmin) : exact_crit_1D(xi){
    this->hmin = hmin;
    this->hmin2 = hmin * hmin;
    outer_diff_abs_1D();
  }

public:
  Eigen::VectorXd compute(Eigen::ArrayXd H){
    Eigen::ArrayXd H_inv = H.inverse();
    Eigen::ArrayXd R2 = (hmin * H_inv).square();
    Eigen::VectorXd pen =  1.875 * (1 + R2 * (R2 - 6) / 21.0) * H_inv / double(n); // 1.875 = 15 / 8

    Eigen::VectorXd loss(H.size());
    for (int no_h = 0; no_h < H.size(); no_h++){
      Rcpp::checkUserInterrupt();
      double h = H(no_h);

      double r = hmin / h;
      double r2 = r * r;
      double r4 = r2 * r2;
      double r5 = r * r4;

      double c_k1 = - 10 / (3584.0 * h);
      double c_k2_0 = 10 / (112.0 * h);
      double c_k2_1 = -c_k1 / r5;

      double c_p2_0 = 18 * r2;
      double c_p2_1 = -6 * r2;

      double c_p2_2 = 16 * (1 + r * (r - 1) * (5 + r * (r - 4)));
      double c_p2_3 = -5 * (1 + r) * (5 + r * (5 * r - 14));
      double c_p2_4 = 3 * (1 + r * (10 + r));
      double c_p2_5 = 5 * (1 + r);

      int u_ind = 0;
      double itu = u(u_ind);

      while ((itu <  h - hmin) & (u_ind < u.size())){// itu is in [0; h - hmin[
        // (Kh * Kh) (itu)
        double z = itu / h;
        double z2 = z * z;
        double z3 = z * z2;
        double z4 = z * z3;

        double p_1 = std::pow(z - 2, 5);
        double p_2 = 16 + 40 * z + 36 * z2 + 10 * z3 + z4;

        double k1 = c_k1 * p_1 * p_2;

        // (Kh * Khmin) (itu)
        p_1 = r4 + 21 * (z2 - 1) * (z2 - 1);
        p_2 = c_p2_0 * z2 + c_p2_1;

        double k2 = c_k2_0 * (p_1 + p_2);

        loss(no_h) += k1 - 2 * k2;

        u_ind++;
        itu = u(u_ind);
      }
      while ((itu <  h + hmin) & (u_ind < u.size())) {// itu is in [h - hmin; h + hmin[
        // (Kh * Kh) (itu)
        double z = itu / h;
        double z2 = z * z;
        double z3 = z * z2;
        double z4 = z * z3;

        double p_1 = std::pow(z - 2, 5);
        double p_2 = 16 + 40 * z + 36 * z2 + 10 * z3 + z4;

        double k1 = c_k1 * p_1 * p_2;

        // (Kh * Khmin) (itu)
        p_1 = std::pow(1 + r - z, 5);
        p_2 = c_p2_2 + c_p2_3 * z + c_p2_4 * z2 + c_p2_5 * z3 + z4;

        double k2 = c_k2_1 * p_1 * p_2;

        loss(no_h) += k1 - 2 * k2;

        u_ind++;
        itu = u(u_ind);
      }
      while ((itu <  2 * h) & (u_ind < u.size())){// itu is in [h + hmin; 2 * h[
        // (Kh * Kh) (itu)
        double z = itu / h;
        double z2 = z * z;
        double z3 = z * z2;
        double z4 = z * z3;

        double p_1 = std::pow(z - 2, 5);
        double p_2 = 16 + 40 * z + 36 * z2 + 10 * z3 + z4;

        double k1 = c_k1 * p_1 * p_2;

        // (Kh * Khmin) (itu) = 0
        loss(no_h) += k1;
        u_ind++;
        itu = u(u_ind);
      }
      // adds the terms s.t. i = j
      double k1_0 = 5.0 / (7.0 * h);
      double k2_0 = c_k2_0 * (r4 + 21 - 6 * r2 );

      loss(no_h) += n * (k1_0 - k2_0);
      loss(no_h) /= n2;
    }

    Eigen::VectorXd crit = loss + pen;
    return(crit);
  }


};














// [[Rcpp::export]]
Eigen::VectorXd crit_GK_exact_1D(Eigen::VectorXd xi, Eigen::ArrayXd H){
  double hmin = H.minCoeff();
  GK_exact_crit_1D crit = GK_exact_crit_1D(xi, hmin);
  
  
  Eigen::VectorXd crit_comp = crit.compute(H);
  
  return(crit_comp);
}





// [[Rcpp::export]]
Eigen::VectorXd crit_EK_exact_1D(Eigen::VectorXd xi, Eigen::ArrayXd H){
  double hmin = H.minCoeff();
  EK_exact_crit_1D crit = EK_exact_crit_1D(xi, hmin);

  
  Eigen::VectorXd crit_comp = crit.compute(H);
  
  return(crit_comp);
}


// [[Rcpp::export]]
Eigen::VectorXd crit_BK_exact_1D(Eigen::VectorXd xi, Eigen::ArrayXd H){
  double hmin = H.minCoeff();
  BK_exact_crit_1D crit = BK_exact_crit_1D(xi, hmin);

  
  Eigen::VectorXd crit_comp = crit.compute(H);
  
  return(crit_comp);
}


class criterion_mD{
public:
  Eigen::MatrixXd xi;
  Eigen::MatrixXd hmin;
  Eigen::MatrixXd hmin2;
  Eigen::VectorXd hmin_diag;
  Eigen::MatrixXd P;// P sert pour la dichotomie
  Eigen::MatrixXd Pinv;
  
protected :
  int d;
  int n;
  double n2;
  
public :
  criterion_mD(Eigen::MatrixXd xi) {
    this->xi = xi;
    this->d = xi.cols();
    this->n = xi.rows();
    this->n2 = n * n;
  }
  
public:
  virtual Eigen::VectorXd compute(List H) = 0;
  
};




class exact_crit_mD : public criterion_mD {
  
protected:
  // u contains all (oi - oj)^2 (diagonal case) or all (oi - oj) (full case ) 
  // for all i and j s.t. j < i
  Eigen::MatrixXd u;
  int u_rows;
  
public :
  exact_crit_mD(Eigen::MatrixXd xi):criterion_mD(xi){
    this->u_rows = n * (n - 1) / 2;
    this->u = Eigen::MatrixXd::Zero(u_rows, d);
  }
  
protected :
  void outer_diff_square_mD(){
    // Computation of (oi - oj)^2 for all i and j s.t. j < i
    for (int no_d = 0; no_d < d; no_d++){
      Eigen::VectorXd xi_col = xi.col(no_d);
      int pos_u = 0;
      int n_u = n - 1;
      for (int i = 0; i < (n - 1); i++){
        Rcpp::checkUserInterrupt();
        u.block(pos_u, no_d, n_u, 1) = (xi_col.segment(i + 1, n_u).array() - xi_col(i)).square();
        pos_u += n_u;
        n_u -= 1;
      }
    }
  }
  
  
  void outer_diff_mD(){
    // Computation of (oi - oj) for all i and j s.t. j < i
    for (int no_d = 0; no_d < d; no_d++){
      Eigen::VectorXd xi_col = xi.col(no_d);
      int pos_u = 0;
      int n_u = n - 1;
      for (int i = 0; i < (n - 1); i++){
        Rcpp::checkUserInterrupt();
        u.block(pos_u, no_d, n_u, 1) = xi_col.segment(i + 1, n_u).array() - xi_col(i);
        pos_u += n_u;
        n_u -= 1;
      }
    }
  }
  
  
  
};





class GK_exact_crit_mD_diag : public exact_crit_mD{
private :
  double dlog2pi;
  double dlog2;
  double c_pen;
  
  
public :
  GK_exact_crit_mD_diag(Eigen::MatrixXd xi) : exact_crit_mD(xi){
    double cst_diag = 1 / (M_SQRT2 * M_SQRT_PI * std::pow(n, 1 / double(d)));
    this->hmin = Eigen::VectorXd::Constant(d, 1, cst_diag);
    this->hmin2 = hmin.array().square();
    outer_diff_square_mD();
    this->dlog2pi = d * M_LN_2PI;
    this->dlog2 = d * M_LN2;
    this->c_pen = std::pow(M_1_SQRT_2PI, d) * 2.0 / double(n);
  }
  
  
  GK_exact_crit_mD_diag(Eigen::MatrixXd xi, Eigen::VectorXd hmin) : exact_crit_mD(xi){
    // double cst_diag = 1 / (M_SQRT2 * M_SQRT_PI * std::pow(n, 1 / double(d)));
    this->hmin = hmin;
    this->hmin2 = hmin.array().square();
    outer_diff_square_mD();
    this->dlog2pi = d * M_LN_2PI;
    this->dlog2 = d * M_LN2;
    this->c_pen = std::pow(M_1_SQRT_2PI, d) * 2.0 / double(n);
  }
  
  
  
  
public:
  Eigen::VectorXd compute(List H){
    int nh = H.size();
    
    Eigen::VectorXd pen(nh);
    Eigen::VectorXd loss(nh);
    
    Eigen::VectorXd den(u_rows);
    Eigen::VectorXd sumCarre(u_rows);
    
    for (int no_h = 0; no_h < nh; no_h++){
      
      Eigen::VectorXd h = H(no_h);
      Eigen::VectorXd h2 = h.array().square();
      Eigen::VectorXd D = h2 + hmin2;
      double det = D.prod();
      pen(no_h) = det;
      
      double den_0 = dlog2 + dlog2pi + 2 * h.array().log().sum();
      
      Eigen::VectorXd inv_D1 = 0.5 * h2.cwiseInverse();
      sumCarre = u * inv_D1;
      den = -0.5 * (sumCarre.array() + den_0);
      den = den.array().exp();
      Eigen::VectorXd K1 = den;
      
      den_0 *= -0.5;
      den_0 = std::exp(den_0);
      double K1_0 = den_0;
      
      den_0 = D.array().log().sum() + dlog2pi;
      Eigen::VectorXd  inv_D2 = D.cwiseInverse();
      
      sumCarre = u * inv_D2;
      
      den = -0.5 * (sumCarre.array() + den_0);
      den = den.array().exp();
      Eigen::VectorXd K2 = den;
      
      den_0 *= -0.5;
      den_0 = std::exp(den_0);
      double K2_0 = den_0;
      
      Eigen::VectorXd K = K1 - 2 * K2;
      double sum_K = K.sum();
      loss(no_h) = 2.0 * sum_K;
      // ajout des termes tq i=j
      loss(no_h) += n * (K1_0 - 2 * K2_0);
    }
    
    pen = pen.cwiseSqrt();
    pen = pen.cwiseInverse() * c_pen;
    
    Eigen::VectorXd crit = pen + loss / n2;
    return(crit);
  }
  
  
};


class GK_exact_crit_mD_full : public exact_crit_mD{
  
  
  
public :
  // GK_exact_crit_mD_full(Eigen::MatrixXd xi, Eigen::MatrixXd S) : exact_crit_mD(xi){
  //   // S is the covariance matrix of xi
  //   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(S);
  //   // this->P = eigensolver.eigenvectors();
  //   // this->Pinv = P.inverse();
  //   Eigen::MatrixXd P = eigensolver.eigenvectors();
  //   Eigen::MatrixXd Pinv = P.inverse();
  //   double cst_diag = 1 / (M_SQRT2 * M_SQRT_PI * std::pow(n, 1 / double(d)));
  //   
  //   this->hmin_diag = Eigen::VectorXd::Constant(d, 1, cst_diag);
  //   Eigen::MatrixXd D = hmin_diag.asDiagonal();//cst_diag * Eigen::MatrixXd::Identity(d, d);
  //   this->hmin = P * D * Pinv;
  //   this->hmin2 = this->hmin * this->hmin;
  //   outer_diff_mD();
  // }
  GK_exact_crit_mD_full(Eigen::MatrixXd xi) : exact_crit_mD(xi){
    // S is the covariance matrix of xi
    Eigen::MatrixXd centered = xi.rowwise() - xi.colwise().mean();
    Eigen::MatrixXd S = (centered.transpose() * centered) / double(n - 1);
    // diagonalisation de S
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(S);
    this->P = eigensolver.eigenvectors();
    this->Pinv = P.inverse();
    // Eigen::MatrixXd P = eigensolver.eigenvectors();
    // Eigen::MatrixXd Pinv = P.inverse();
    double cst_diag = 1 / (M_SQRT2 * M_SQRT_PI * std::pow(n, 1 / double(d)));
    
    this->hmin_diag = Eigen::VectorXd::Constant(d, 1, cst_diag);
    Eigen::MatrixXd D = hmin_diag.asDiagonal();//cst_diag * Eigen::MatrixXd::Identity(d, d);
    this->hmin = P * D * Pinv;
    this->hmin2 = this->hmin * this->hmin;
    outer_diff_mD();
  }
  
  
  // GK_exact_crit_mD_full(Eigen::MatrixXd xi, Eigen::MatrixXd S, Eigen::VectorXd hmin_diag) : exact_crit_mD(xi){
  GK_exact_crit_mD_full(Eigen::MatrixXd xi, Eigen::MatrixXd hmin) : exact_crit_mD(xi){
  // GK_exact_crit_mD_full(Eigen::MatrixXd xi, Eigen::MatrixXd S, List H_init) : exact_crit_mD(xi){
    // S is the covariance matrix of xi
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(S);
    // this->P = eigensolver.eigenvectors();
    // this->Pinv = P.inverse();
    // // double cst_diag = 1 / (M_SQRT2 * M_SQRT_PI * std::pow(n, 1 / double(d)));
    // // Eigen::MatrixXd D = cst_diag * Eigen::MatrixXd::Identity(d, d);
    // 
    // // Pour chaque élément de H_init on applique le changement de base lié à la matrice de passage P
    // List H;
    // for (int i=0; i<H_init.size(); i++){
    //   Eigen::VectorXd diag_H = H_init(i);
    //   Eigen::MatrixXd Hi_tmp = P * diag_H.asDiagonal() * Pinv;
    //   H.push_back(Hi_tmp);
    // }
    // this->hmin_diag = Eigen::VectorXd::Constant(d, 1, terme_diag_hmin);
    // Eigen::MatrixXd D = hmin_diag.asDiagonal();
    // this->hmin = P * D * Pinv;
    this->hmin = hmin;
    this->hmin2 = hmin * hmin;
    
    outer_diff_mD();
  }
  
  
  
public:
  Eigen::VectorXd compute(List H){
    int nh = H.size();
    // Rcout << "crit compute" << std::endl;
    // penalty computation
    double cst = 2.0 * std::pow(M_1_SQRT_2PI, d) / double(n);
    Eigen::VectorXd pen = Eigen::VectorXd::Constant(nh, cst);
    
    // loss computation
    Eigen::VectorXd loss = Eigen::VectorXd::Zero(nh);
    double s2 = d * M_LN_2PI;
    for (int no_h = 0; no_h < nh; no_h++){
      Rcpp::checkUserInterrupt();
      // Rcout << "no_h = " << no_h << std::endl;
      Eigen::MatrixXd h = H(no_h);
      Eigen::MatrixXd h2 = h * h;
      
      
      // Computation of the multivariate gaussian density K_{sqrt(2)h}
      // For a faster computation the cholesky decomposition of sig is used
      
      Eigen::MatrixXd sig2 = 2 * h2;
      // the matrix of the decomposition
      Eigen::MatrixXd L_1 = sig2.llt().matrixL();
      // solve the linear matricial equation Lx = u
      Eigen::MatrixXd x_1 = L_1.colPivHouseholderQr().solve(u.transpose());
      
      Eigen::VectorXd LL_1 = L_1.diagonal().array().log();
      double s1_1 = LL_1.sum(); // = log(det(sig))
      
      Eigen::MatrixXd sum_1 = x_1.array().square();
      Eigen::VectorXd den_1 = sum_1.colwise().sum();
      
      den_1 = den_1.array() + 2 * s1_1 + s2;
      den_1 = - den_1 / 2.0;
      den_1 = den_1.array().exp();
      Eigen::VectorXd K1 = den_1;
      double K1_0 = std::exp(-0.5 * (s2 + 2 * s1_1));
      
      // Computation of the multivariate gaussian density K_{sqrt(h^2 + hmin^2)}
      sig2 = h2 + hmin2;
      double det = sig2.determinant();
      
      pen(no_h) /= std::sqrt(det);
      
      Eigen::MatrixXd L_2 = sig2.llt().matrixL();
      // solve the linear matricial equation Lx = u
      Eigen::MatrixXd x_2 = L_2.colPivHouseholderQr().solve(u.transpose());
      
      Eigen::VectorXd LL_2 = L_2.diagonal().array().log();
      double s1_2 = LL_2.sum();
      
      Eigen::MatrixXd sum_2 = x_2.array().square();
      Eigen::VectorXd den_2 = sum_2.colwise().sum();
      
      den_2 = den_2.array() + 2 * s1_2 + s2;
      den_2 = - den_2 / 2.0;
      den_2 = den_2.array().exp();
      Eigen::VectorXd K2 = den_2;
      double K2_0 = std::exp(-0.5 * (s2 + 2 * s1_2));
      
      Eigen::VectorXd K = K1 - 2 * K2;
      double K0 = K1_0 - 2 * K2_0;
      
      loss(no_h) = 2 * K.sum();
      
      // add terme such that i=j
      loss(no_h) += n * K0;
    }
    loss /= n2;
    Eigen::VectorXd crit = loss + pen;
    return(crit);
  }
  
  
};













//crit_GK_exact_mD_diag_17
// [[Rcpp::export]]
Eigen::VectorXd crit_GK_exact_mD_diag(Eigen::MatrixXd x_i, List H, Eigen::VectorXd hmin){
  GK_exact_crit_mD_diag crit = GK_exact_crit_mD_diag(x_i, hmin);
  
  // Rcout << "crit.u = " << crit.u << std::endl;
  Eigen::VectorXd crit_comp = crit.compute(H);
  // Eigen::VectorXd crit_comp = Eigen::VectorXd::Zero(H.size());
  return(crit_comp);
}



//crit_GK_exact_mD_full_18
// [[Rcpp::export]]
Eigen::VectorXd crit_GK_exact_mD_full(Eigen::MatrixXd x_i, List H, Eigen::MatrixXd hmin){
  GK_exact_crit_mD_full crit = GK_exact_crit_mD_full(x_i, hmin);
  
  
  Eigen::VectorXd crit_comp = crit.compute(H);
  // Rcout << "crit.u = " << crit.u << std::endl;
  // Eigen::VectorXd crit_comp = crit.compute(H);
  // Eigen::VectorXd crit_comp = Eigen::VectorXd::Zero(H.size());
  // List l = dich_mD_diag(crit, nh_max, tol);
  
  // Eigen::VectorXd h_opt = dich_mD_diag(crit, nh_max, tol);
  // Eigen::MatrixXd h_opt = Eigen::MatrixXd::Zero(2,2);
  return(crit_comp);
}

























// risque Lp avec le noyau gaussien en 1D
// [[Rcpp::export]]
double risque_Lp(int p, NumericVector f, NumericVector o_i, 
                        NumericVector x, double hopt_c, double a, double b){
  
  
  //int nc = 1;
  int nQMC = x.size();
  int n = o_i.size();
  
  double risk;
  
  const double sqrt_2Pi = std::sqrt(2.0 * PI);
  double n_sqrt_2Pi = double(n)*sqrt_2Pi;
  
  
  
    double h = hopt_c;
    double g = 0.0;
    
    for (NumericVector::iterator itx = x.begin(), itf = f.begin();
         itx != x.end(); ++itx, ++itf){
      
      double xk = *itx;
      double somme_i = 0.0;
      
      for (NumericVector::iterator ito = o_i.begin(); ito != o_i.end(); ++ito){
        //(cas du noyau gaussien) les autres noyaux sont à faire
        double oi = *ito;
        double u = (oi - xk) / h;
        double Kh = std::exp(-u * u / 2.0); 
        somme_i += Kh; // somme sur i
      }
      
      double hat_f_i = somme_i / (h * n_sqrt_2Pi);
      double fk = *itf;
      g += std::pow(std::abs(hat_f_i - fk), p); // somme sur k de la val absolue à la puissance p
    }
    g *= (b - a) / double(nQMC) ;
    
    risk = g;
  // }
  
  
  return(risk);
  
}















// risque Lp avec le noyau gaussien en multidim
// la boucle sur n est remplacee par un calcul vectoriel pour le calcul de la somme des Kh(x-o_i)
// vectorisation de la boucle sur no_d
// [[Rcpp::export]]
double risque_Lp_mD_3(int p, Eigen::VectorXd f, Eigen::MatrixXd o_i, int nQMC, 
                      Eigen::MatrixXd h, Eigen::VectorXd a, Eigen::VectorXd b){
  
  //int nc = 1;
  
  int nloc = nQMC;// = nloc pour avoir tous les points il faut faire un outer de toutes les composantes 
  int n = o_i.rows();
  int d = o_i.cols();
  
  double risk=0;
  
  
  
  
  
  
  
  //const double sqrt_2Pi = M_SQRT2 * M_SQRT_PI;//std::sqrt(2.0 * PI); 
  //double n_sqrt_2Pi = double(n) * sqrt_2Pi;
  double dlog2Pi = d * M_LN_2PI;
  double dlog2PiD2 = 0.5 * dlog2Pi;
  
  
  
  Eigen::MatrixXd h2 = h * h;
  
  double g = 0;
  
  int size = std::pow(nloc, d);
  Eigen::ArrayXd point = Eigen::VectorXd::Zero(d);
  Eigen::VectorXd xi(d);
  
  Eigen::MatrixXd u(n, d);
  
  Eigen::ArrayXd delta = (b - a) / (nloc - 1);
  
  
  Eigen::ArrayXi pows_d = Eigen::ArrayXi::LinSpaced(d, 0, d - 1);
  Eigen::ArrayXi denominators = pow(nloc, pows_d);
  
  Eigen::ArrayXi quot(d);
  
  
  // the matrix of the decomposition : h2 = L*t(L) 
  Eigen::MatrixXd L_1 = h2.llt().matrixL();
  
  
  Eigen::VectorXd LL_1 = L_1.diagonal().array().log();
  double log_det_h = LL_1.sum(); // = log(det(h)) = 0.5 * log(det(sig))
  
  
  double c2 = - log_det_h - dlog2PiD2;
    
    
  Eigen::MatrixXd x_1(n, d), squared(n, d);
  
  Eigen::VectorXd txx(n), den(n);
  
  
  double hat_f, fx, somme_i;
  
  
  
  for (int no_x=0; no_x<size; no_x++){
    Rcpp::checkUserInterrupt();
   
    quot = no_x / denominators;// fait une division entiere 
    
    point = (quot - (quot / nloc) * nloc).cast<double>();
    
    xi = a + (point * delta).matrix();
    
    
    
     
    u = o_i.rowwise() - xi.transpose();
    
    // solve the linear matricial equation Lx = u
    x_1 = L_1.colPivHouseholderQr().solve(u.transpose());
    
    // calcul de la diagonale de t(x_1) * x_1
    
    txx = x_1.array().square().colwise().sum();
    
    den = - 0.5 * txx.array() + c2;
    
    somme_i = den.array().exp().sum();
    
    hat_f = somme_i / double(n);// / (h * n_sqrt_2Pi);
    
    fx = f(no_x);
    g += std::pow(std::abs(hat_f - fx), p); // somme sur k de la val absolue à la puissance p
  }
  
  
  
  g *= (b - a).prod() / double(size) ;
  
  risk = g;
  
  
  
  
  
  
  return(risk);
  
}
