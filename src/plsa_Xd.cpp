// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

  

// Normalize a matrix
void rescale_back2_prob(Eigen::MatrixXd & pM,double threshold){
  // loop over and make sure we don't have zeros
  auto ptr = pM.data();
  for (int i = 0; i < pM.size(); i++){
    *(ptr+i) = (*(ptr+i) >= threshold) ? *(ptr+i) : threshold;
  }
  
  // re-scale
  Eigen::VectorXd mSum = pM.colwise().sum().cwiseInverse();
  pM = pM*mSum.asDiagonal();

   
}


void rescale_back2_prob(Eigen::VectorXd & pV,double threshold){
  // loop over and make sure we don't have zeros
  auto ptr = pV.data();
  for (int i = 0; i < pV.size(); i++){
    *(ptr+i) = (*(ptr+i) >= threshold) ? *(ptr+i) : threshold;
  }
  
  double inv_sum =1.0/ pV.sum();
  pV *= inv_sum;
}


// [[Rcpp::export]]
Rcpp::List sparse_plsa_3d(const Eigen::VectorXi & aye,
                   const Eigen::VectorXi & jay,
                   const Eigen::VectorXi & kay,
                   const Eigen::VectorXd & qty,
                   const int num_a,
                   const int num_b,
                   const int num_c,
                   const int nz,
                   const int niter){

  // Initialize Matrices
  auto num = qty.size();
  double threshhold = 0.1/double(num);

  // Initialize vectors
  Eigen::VectorXd conv(niter);
  Eigen::VectorXd Pz = Eigen::VectorXd::Zero(nz);
  Eigen::VectorXd tPz(nz);
  Eigen::VectorXd P_abc(num);
  Eigen::MatrixXd Pa_z = Eigen::MatrixXd::Random(num_a,nz);
  Eigen::MatrixXd Pb_z = Eigen::MatrixXd::Random(num_b,nz);
  Eigen::MatrixXd Pc_z = Eigen::MatrixXd::Random(num_c,nz);
  
  Eigen::MatrixXd tPa_z(num_a,nz); // temporary storage
  Eigen::MatrixXd tPb_z(num_b,nz);
  Eigen::MatrixXd tPc_z(num_c,nz);

/*
  if(initialization == 1){
    Eigen::VectorXi Pi_a = Eigen::VectorXi::Random(num_a);
    Eigen::VectorXi Pi_b = Eigen::VectorXi::Random(num_b);
    Eigen::VectorXi Pi_c = Eigen::VectorXi::Random(num_c);
    Pi_a = Pi_a.cwiseAbs();
    Pi_b = Pi_b.cwiseAbs();
    Pi_c = Pi_c.cwiseAbs();
    
    for(auto i=0;i<num;i++){
      int a = aye(i)-1;
      int b = jay(i)-1;
      int c = kay(i)-1;
      int za = Pi_a(a)%nz;
      int zb = Pi_b(b)%nz;
      int zc = Pi_c(c)%nz;
      double temp = qty(i);
      Pa_z(a,za) += temp;
      Pb_z(b,zb) += temp;
      Pc_z(c,zc) += temp;
    }

  }
  else{
 */
  Pa_z = Pa_z.array().abs();
  Pb_z = Pb_z.array().abs();
  Pc_z = Pc_z.array().abs();


  rescale_back2_prob(Pa_z,0.1);
  rescale_back2_prob(Pb_z,0.1);
  rescale_back2_prob(Pc_z,0.1);
  rescale_back2_prob(Pz,1.);

  //===================================
  // Iterate the EM algorithm
  //===================================
  for(int it=0;it<niter;it++){
  
    //------------------------
    // E step
    //------------------------
    double cost=0;
    // loop over data
    for(auto i=0;i<num;i++){
      int a = aye(i)-1;// R uses 1-based indices
      int b = jay(i)-1;
      int c = kay(i)-1;

      // loop over latent variables
      double tP_abc = 0;
      for(int z=0;z<nz;z++){
        tP_abc += Pz(z)*Pa_z(a,z)*Pb_z(b,z)*Pc_z(c,z);
      }
      
      // don't divide by zero
      tP_abc = (tP_abc<threshhold)?threshhold:tP_abc;
      
      // Update cost
      cost -= log(tP_abc)*qty(i);
      
      // store updated tensor fit
      P_abc(i) = tP_abc;
    }
    conv[it] = cost;

    
    //------------------------
    // M step
    //------------------------
    tPa_z.setZero();
    tPb_z.setZero();
    tPc_z.setZero();
    tPz.setZero();
    
    // loop over data,    
    for(auto i=0;i<num;i++){
      int a = aye(i)-1;
      int b = jay(i)-1;
      int c = kay(i)-1;
      double commonRatio = qty(i)/P_abc(i);
      for(int z=0;z<nz;z++){
        double temp = commonRatio*Pz(z)*Pa_z(a,z)*Pb_z(b,z)*Pc_z(c,z);
        tPa_z(a,z) += temp;
        tPb_z(b,z) += temp;
        tPc_z(c,z) += temp;
        tPz(z) += temp;
      }
    }
    
    //------------------------
    //Copy and normalize
    //------------------------
    Pa_z = tPa_z;
    Pb_z = tPb_z;
    Pc_z = tPc_z;
    Pz = tPz;
    
    rescale_back2_prob(Pa_z,threshhold);
    rescale_back2_prob(Pb_z,threshhold);
    rescale_back2_prob(Pc_z,threshhold);
    rescale_back2_prob(Pz,threshhold);
    

  }

  return Rcpp::List::create(Rcpp::Named("Pa_z")=Pa_z,
                            Rcpp::Named("Pb_z")=Pb_z,
                            Rcpp::Named("Pc_z")=Pc_z,
                            Rcpp::Named("Pz")=Pz,
                            Rcpp::Named("cost")=conv);

}
  

 
 


// PLSA 

// [[Rcpp::export]]
Rcpp::List sparse_plsa_2d(const Eigen::VectorXi & aye,
                          const Eigen::VectorXi & jay,
                          const Eigen::VectorXd & qty,
                          const int num_a,
                          const int num_b,
                          const int nz,
                          const int niter){
  
  
  
  // Get sizes
  
  // Initialize Matrices
  auto num = qty.size();
  double threshhold = 0.1/double(num);
  Eigen::VectorXd conv(niter);
  
  // Initialize vectors
  Eigen::VectorXd Pz = Eigen::VectorXd::Zero(nz);
  Eigen::VectorXd tPz(nz);
  Eigen::VectorXd P_ab(num);
  Eigen::MatrixXd Pa_z = Eigen::MatrixXd::Random(num_a,nz);
  Eigen::MatrixXd Pb_z = Eigen::MatrixXd::Random(num_b,nz);

  Eigen::MatrixXd tPa_z(num_a,nz);
  Eigen::MatrixXd tPb_z(num_b,nz);

  
  //normalize
  Pa_z = Pa_z.array().abs();
  Pb_z = Pb_z.array().abs();
  
  
  rescale_back2_prob(Pa_z,0.1);
  rescale_back2_prob(Pb_z,0.1);
  rescale_back2_prob(Pz,1.);
  
  
  //===================================
  // Iterate the EM algorithm
  //===================================
  for(int it=0;it<niter;it++){
    
    //------------------------
    // E step
    //------------------------
    
    double cost=0;
    
    // loop over data
    for(auto i=0;i<num;i++){
      int a = aye(i)-1;// R uses 1-based indices
      int b = jay(i)-1;
      double tP_ab = 0;
      
      
      
      // loop over latent variables
      for(int z=0;z<nz;z++){
        tP_ab += Pz(z)*Pa_z(a,z)*Pb_z(b,z);
        
      }
      // don't divide by zero
      tP_ab = (tP_ab<threshhold)?threshhold:tP_ab;
      
      // Update cost
      cost -= log(tP_ab)*qty(i);
      P_ab(i) = tP_ab;
        
    }
    conv[it] = cost;
    
    
    //------------------------
    // M step
    //------------------------
    tPa_z.setZero();
    tPb_z.setZero();
    tPz.setZero();
    
    // loop over data,    
    for(auto i=0;i<num;i++){
      int a = aye(i)-1;
      int b = jay(i)-1;
      double commonRatio = qty(i)/P_ab(i);
      for(int z=0;z<nz;z++){
        double temp = commonRatio*Pz(z)*Pa_z(a,z)*Pb_z(b,z);
        tPa_z(a,z) += temp;
        tPb_z(b,z) += temp;
        tPz(z) += temp;
      }
    }
    
    //------------------------
    //Copy and normalize
    //------------------------
    Pa_z = tPa_z;
    Pb_z = tPb_z;
    Pz = tPz;
    
    rescale_back2_prob(Pa_z,threshhold);
    rescale_back2_prob(Pb_z,threshhold);
    rescale_back2_prob(Pz,threshhold);
    
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Pa_z")=Pa_z,
                            Rcpp::Named("Pb_z")=Pb_z,
                            Rcpp::Named("Pz")=Pz,
                            Rcpp::Named("cost")=conv);
  
}

//======================================
// 4 Way PLSA
//======================================

// [[Rcpp::export]]
Rcpp::List sparse_plsa_4d(const Eigen::VectorXi & aye,
                          const Eigen::VectorXi & jay,
                          const Eigen::VectorXi & kay,
                          const Eigen::VectorXi & ell,
                          const Eigen::VectorXd & qty,
                          const int num_a,
                          const int num_b,
                          const int num_c,
                          const int num_d,
                          const int nz,
                          const int niter){
  
  // Initialize Matrices
  auto num = qty.size();
  double threshhold = 0.1/double(num); // don't allow numbers to fall below this
  
  // Initialize vectors
  Eigen::VectorXd conv(niter);
  Eigen::VectorXd Pz = Eigen::VectorXd::Zero(nz);
  Eigen::VectorXd tPz(nz);
  Eigen::VectorXd P_abcd(num);
  Eigen::MatrixXd Pa_z = Eigen::MatrixXd::Random(num_a,nz);
  Eigen::MatrixXd Pb_z = Eigen::MatrixXd::Random(num_b,nz);
  Eigen::MatrixXd Pc_z = Eigen::MatrixXd::Random(num_c,nz);
  Eigen::MatrixXd Pd_z = Eigen::MatrixXd::Random(num_d,nz);
  
  Eigen::MatrixXd tPa_z(num_a,nz); // temporary storage
  Eigen::MatrixXd tPb_z(num_b,nz);
  Eigen::MatrixXd tPc_z(num_c,nz);
  Eigen::MatrixXd tPd_z(num_d,nz);
  
  // Remove negatives and normalize
  Pa_z = Pa_z.array().abs();
  Pb_z = Pb_z.array().abs();
  Pc_z = Pc_z.array().abs();
  Pd_z = Pd_z.array().abs();
  rescale_back2_prob(Pa_z,0.1);
  rescale_back2_prob(Pb_z,0.1);
  rescale_back2_prob(Pc_z,0.1);
  rescale_back2_prob(Pd_z,0.1);
  rescale_back2_prob(Pz,1.);
  
  //===================================
  // Iterate the EM algorithm
  //===================================
  for(int it=0;it<niter;it++){
    
    //------------------------
    // E step
    //------------------------
    double cost=0;
    // loop over data
    for(auto i=0;i<num;i++){
      int a = aye(i)-1;// R uses 1-based indices
      int b = jay(i)-1;
      int c = kay(i)-1;
      int d = ell(i)-1;
      
      // loop over latent variables
      double tP_abcd = 0;
      for(int z=0;z<nz;z++){
        tP_abcd += Pz(z)*Pa_z(a,z)*Pb_z(b,z)*Pc_z(c,z)*Pd_z(d,z);
      }
      
      // don't divide by zero
      tP_abcd = (tP_abcd<threshhold)?threshhold:tP_abcd;
      
      // Update cost
      cost -= log(tP_abcd)*qty(i);
      
      // store updated tensor fit
      P_abcd(i) = tP_abcd;
    }
    conv[it] = cost;
    
    
    //------------------------
    // M step
    //------------------------
    tPa_z.setZero();
    tPb_z.setZero();
    tPc_z.setZero();
    tPd_z.setZero();
    tPz.setZero();
    
    // loop over data,    
    for(auto i=0;i<num;i++){
      int a = aye(i)-1;
      int b = jay(i)-1;
      int c = kay(i)-1;
      int d = ell(i)-1;
      
      double commonRatio = qty(i)/P_abcd(i);
      for(int z=0;z<nz;z++){
        double temp = commonRatio*Pz(z)*Pa_z(a,z)*Pb_z(b,z)*Pc_z(c,z)*Pd_z(d,z);
        tPa_z(a,z) += temp;
        tPb_z(b,z) += temp;
        tPc_z(c,z) += temp;
        tPd_z(d,z) += temp;
        tPz(z) += temp;
      }
    }
    
    //------------------------
    //Copy and normalize
    //------------------------
    Pa_z = tPa_z;
    Pb_z = tPb_z;
    Pc_z = tPc_z;
    Pd_z = tPd_z;
    Pz = tPz;
    
    rescale_back2_prob(Pa_z,threshhold);
    rescale_back2_prob(Pb_z,threshhold);
    rescale_back2_prob(Pc_z,threshhold);
    rescale_back2_prob(Pd_z,threshhold);
    rescale_back2_prob(Pz,threshhold);
    
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Pa_z")=Pa_z,
                            Rcpp::Named("Pb_z")=Pb_z,
                            Rcpp::Named("Pc_z")=Pc_z,
                            Rcpp::Named("Pd_z")=Pd_z,
                            Rcpp::Named("Pz")=Pz,
                            Rcpp::Named("cost")=conv);
  
}


  
  
  
