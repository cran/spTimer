
#include "main_gpp.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"
//#include "Print.h"



// with zfit summary values (mean and variance/sd)
// with one phi parameter
void GIBBS_zfitsum_onephi_gpp(int *cov, int *spdecay, double *flag, int *its, int *burnin,
     int *n, int *m, int *T, int *r, int *rT, int *p, int *N, int *report,
	 int *ft,
     double *shape_e, double *shape_eta, double *shape_l,  
     double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *mu_beta, double *delta2_beta,
     double *mu_rho,  double *delta2_rho, double *alpha_l, double *delta2_l,
     double *phi_eta, double *tau_eta, double *phis, int *phik, 
     double *dm, double *dnm, int *constant, 
     double *sig2e, double *sig2eta, double *sig2l, double *beta, 
     double *rho, double *mu_l, double *X, double *z, double *w0, double *w,
     int *transform, double *phi_etaf, double *accept_etaf, double *nupf,
     double *sig2ef, double *sig2etaf, double *betaf, double *rhof, 
     double *mu_lf, double *sig2lf, double *w0f, double *wf, 
     double *gof, double *penalty, double *z_mean_sd)
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, col, i, j, m1, r1, T1, p1, N1, rep1, brin;
     its1 = *its;
     col = *constant;
//     n1 = *n;
     m1 = *m;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     N1 = *N;
     rep1 = *report;
     brin = *burnin;
//     trans1 = *transform;

     double accept1; //, mn_rep[N1], var_rep[N1];
     accept1 = 0.0;
     double *mn_rep, *var_rep;
     mn_rep = (double *) malloc((size_t)((N1)*sizeof(double)));          
     var_rep = (double *) malloc((size_t)((N1)*sizeof(double)));          

     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      

     double *phi_etap, *sig2ep, *sig2etap, *rhop, *betap;
     double *mu_lp, *sig2lp, *wp, *w0p, *zfit;
     double *phi_eta1, *sig2e1, *sig2eta1, *rho1, *beta1;
     double *mu_l1, *sig2l1, *w01, *w1;
     double *acc_eta, *out; // *z1;
          
     phi_etap = (double *) malloc((size_t)((1)*sizeof(double)));          
     sig2ep = (double *) malloc((size_t)((1)*sizeof(double)));
     sig2etap = (double *) malloc((size_t)((1)*sizeof(double)));
     rhop = (double *) malloc((size_t)((1)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double)));          
     mu_lp = (double *) malloc((size_t)((r1)*sizeof(double)));    
     sig2lp = (double *) malloc((size_t)((r1)*sizeof(double))); 
     wp = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double))); 
     w0p = (double *) malloc((size_t)((m1*r1)*sizeof(double))); 
     zfit = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi_eta1 = (double *) malloc((size_t)((1)*sizeof(double)));
     sig2e1 = (double *) malloc((size_t)((1)*sizeof(double)));
     sig2eta1 = (double *) malloc((size_t)((1)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((1)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));     
     mu_l1 = (double *) malloc((size_t)((r1)*sizeof(double)));    
     sig2l1 = (double *) malloc((size_t)((r1)*sizeof(double))); 
     w01 = (double *) malloc((size_t)((m1*r1)*sizeof(double))); 
     w1 = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double))); 
     
     acc_eta = (double *) malloc((size_t)((1)*sizeof(double)));
     out = (double *) malloc((size_t)((1)*sizeof(double)));     
//     z1 = (double *) malloc((size_t)((N1)*sizeof(double)));

     double *nu, *nup;
     nu = (double *) malloc((size_t)((1)*sizeof(double)));
     nup = (double *) malloc((size_t)((1)*sizeof(double)));     
     nu[0] = 2.0/3.0;
          
       ext_sige(phi_eta, phi_eta1);
       ext_sige(sig2e, sig2e1);
       ext_sigeta(sig2eta, sig2eta1);
       ext_rho(rho, rho1);
       ext_beta(p, beta, beta1);
       ext_mul(r, mu_l, mu_l1);
       ext_sigl(r, sig2l, sig2l1);
//       ext_o(N, z, z1);
       for(j=0; j<m1*r1*T1; j++){
         w1[j] = w[j];
       }     
       for(j=0; j<m1*r1; j++){
         w01[j] = w0[j];
       }     
            
/*
// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          out[0]=z[j];
          mvrnormal(constant, out, sig2e1, constant, out);  
          z1[j] =  out[0];
          }
          else {
          z1[j] = z[j];
          }     
     }
*/

     GetRNGstate();
     for(i=0; i < its1; i++) {

     JOINT_onephi_gpp(cov, spdecay, flag, n, m, T, r, rT, p, N, 
     shape_e, shape_eta, shape_l, 
     phi_a, phi_b,
     prior_a, prior_b, mu_beta, delta2_beta, 
     mu_rho, delta2_rho, alpha_l, delta2_l, phi_eta1, tau_eta, phis, phik, nu, 
     dm, dnm, constant, sig2e1, sig2eta1, sig2l1, beta1, rho1, mu_l1, X, z, 
     w01, w1, phi_etap, acc_eta, nup, sig2ep, sig2etap, betap, rhop, mu_lp, sig2lp, 
     w0p, wp, zfit);

     accept1 += acc_eta[0];
        
     phi_etaf[i] = phi_etap[0];
     nupf[i] = nup[0];
     
     sig2ef[i] = sig2ep[0];
     sig2etaf[i] = sig2etap[0];
     rhof[i] = rhop[0];
     for(j=0; j < r1; j++) {
         sig2lf[j+i*r1] = sig2lp[j];     
     }
     for(j=0; j < p1; j++){
         betaf[j+i*p1] = betap[j];
     }         
     for(j=0; j < r1; j++) {
         mu_lf[j+i*r1] = mu_lp[j];     
     }
     for(j=0; j<m1*r1; j++){
         w0f[j+i*m1*r1] = w0p[j];
     }     
     for(j=0; j<m1*r1*T1; j++){
         wf[j+i*m1*r1*T1] = wp[j];
     }     


// for pmcc, fitted sum values    
     for(j=0; j < N1; j++){
          if(i >= brin){    
          // Three options: ft: 0=NONE, 1=SQRT, 2=LOG
		  if(ft[0]==0){
              mn_rep[j] += zfit[j];
              var_rep[j] += zfit[j]*zfit[j];
		  }
		  else{
			  if(ft[0]==1){
              mn_rep[j] += zfit[j]*zfit[j];
              var_rep[j] += zfit[j]*zfit[j]*zfit[j]*zfit[j];
			  }
			  else{
              mn_rep[j] += exp(zfit[j]);
              var_rep[j] += exp(zfit[j])*exp(zfit[j]);
			  }
		  }
          }
     }

	 
//	 
       ext_sige(phi_etap, phi_eta1);
       ext_sige(nup, nu);       
       ext_sige(sig2ep, sig2e1);
       ext_sigeta(sig2etap, sig2eta1); 
       ext_rho(rhop, rho1); 
       ext_beta(p, betap, beta1);              
       for(j=0; j<m1*r1*T1; j++){
         w1[j] = wp[j];
       }     
       for(j=0; j<m1*r1; j++){
         w01[j] = w0p[j];
       }     

     if(cov[0]==4){
     para_printRnu (i, its1, rep1, p1, accept1, phi_etap, nup, rhop, sig2ep, sig2etap, betap); 
     }                   
     else {              
     para_printR (i, its1, rep1, p1, accept1, phi_etap, rhop, sig2ep, sig2etap, betap); 
     }
//     printR(i, its1); 
     
     } // end of iteration loop
     PutRNGstate();

     free(phi_etap); free(nu); free(nup); free(sig2ep); free(sig2etap); 
     free(rhop); free(betap); free(mu_lp); free(sig2lp); free(wp); free(w0p); 
     free(zfit); free(phi_eta1); free(sig2e1); free(sig2eta1); free(rho1); 
     free(beta1); free(mu_l1); free(sig2l1); free(acc_eta); free(out); 
     //free(z1);  
     free(w01); free(w1); 

     accept_etaf[0] = accept1;

	int *iit;
    iit = (int *) malloc((size_t)((col)*sizeof(int)));
	iit[0] = its[0] - burnin[0];

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit[0];
          var_rep[j] = var_rep[j]/iit[0];
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          z_mean_sd[j] = mn_rep[j];
          z_mean_sd[N1+j] = sqrt(var_rep[j]);
     }

     double pen, go;
     pen = 0;
     go =0;

// pmcc          
     for(j=0; j < N1; j++){
         if (flag[j] == 1.0){
          mn_rep[j] = 0.0;
          var_rep[j] = 0.0;
         }
         else{
          mn_rep[j] = mn_rep[j];
          var_rep[j] = var_rep[j];
          mn_rep[j] = (mn_rep[j] - z[j])*(mn_rep[j] - z[j]);          
         }
         pen += var_rep[j]; 
         go += mn_rep[j];
     }
     penalty[0] = pen;
     gof[0] = go;       

     free(mn_rep); free(var_rep); free(iit);
     
     return;
}



// for spatially varying GPP
// with zfit summary values (mean and variance/sd)
// with one phi parameter
void GIBBSsp_zfitsum_onephi_gpp(int *intercept, int *cov, int *spdecay, 
     double *flag, int *its, int *burnin, int *n, int *m, int *T, int *r, 
     int *rT, int *p, int *q, int *N, int *report, double *shape_e, 
     double *shape_eta, double *shape_beta, double *shape_l, double *prior_a, 
     double *prior_b, double *mu_beta, double *delta2_beta, double *mu_rho,  
     double *delta2_rho, double *alpha_l, double *delta2_l, double *phi_eta, 
     double *tau_eta, double *phis, int *phik, double *dm, double *dnm, 
     int *constant, 
     double *sig2e, double *sig2eta, double *sig2beta, double *sig2l, 
     double *beta, double *betas, double *rho, double *mu_l, double *X, 
     double *Xsp, double *z, double *w0, double *w, int *transform, 
     double *phi_etaf, double *accept_etaf, double *nupf, double *sig2ef, 
     double *sig2etaf, double *sig2betaf, double *betaf, double *betaspf, 
     double *rhof, double *mu_lf, double *sig2lf, double *w0f, double *wf, 
     double *gof, double *penalty, double *z_mean_sd)
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, col, i, j, m1, r1, T1, p1, q1, N1, rep1, brin;
     its1 = *its;
     col = *constant;
//     n1 = *n;
     m1 = *m;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     q1 = *q;
     N1 = *N;
     rep1 = *report;
     brin = *burnin;
//     trans1 = *transform;


     double accept1; 
     accept1 = 0.0;
     double *mn_rep, *var_rep;
     mn_rep = (double *) malloc((size_t)((N1)*sizeof(double)));          
     var_rep = (double *) malloc((size_t)((N1)*sizeof(double)));          

     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      

     double *phi_etap, *sig2ep, *sig2etap, *sig2betap, *rhop, *betap, *betasp;
     double *mu_lp, *sig2lp, *wp, *w0p, *zfit;
     double *phi_eta1, *sig2e1, *sig2eta1, *sig2beta1, *rho1, *beta1;
     double *mu_l1, *sig2l1; //, *w01, *w1;
     double *acc_eta, *out; // *z1;
          
     phi_etap = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig2ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2etap = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2betap = (double *) malloc((size_t)((col)*sizeof(double)));
     rhop = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1*col)*sizeof(double)));          
     betasp = (double *) malloc((size_t)((q1*m1)*sizeof(double)));          
     mu_lp = (double *) malloc((size_t)((r1*col)*sizeof(double)));    
     sig2lp = (double *) malloc((size_t)((r1*col)*sizeof(double))); 
     wp = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double))); 
     w0p = (double *) malloc((size_t)((m1*r1)*sizeof(double))); 
     zfit = (double *) malloc((size_t)((N1*col)*sizeof(double)));

     phi_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2beta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1*col)*sizeof(double)));     
     mu_l1 = (double *) malloc((size_t)((r1*col)*sizeof(double)));    
     sig2l1 = (double *) malloc((size_t)((r1*col)*sizeof(double))); 
     
     acc_eta = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));     

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 2.0/3.0;
          
       ext_sige(phi_eta, phi_eta1);
       ext_sige(sig2e, sig2e1);
       ext_sigeta(sig2eta, sig2eta1);
       ext_sigeta(sig2beta, sig2beta1);
       ext_rho(rho, rho1);
       ext_beta(p, beta, beta1);
       ext_mul(r, mu_l, mu_l1);
       ext_sigl(r, sig2l, sig2l1);

     GetRNGstate();
     for(i=0; i < its1; i++) {


     JOINT_onephi_sp_gpp(intercept, cov, spdecay, flag, n, m, T, r, rT, p, q,
     N, shape_e, shape_eta, shape_beta, shape_l, prior_a, prior_b, mu_beta, 
     delta2_beta, mu_rho, delta2_rho, alpha_l, delta2_l, phi_eta1, tau_eta,
     phis, phik, nu, dm, dnm, constant, sig2e1, sig2eta1, sig2beta1, sig2l1,
     beta1, betas, rho1, mu_l1, X, Xsp, z, w0, w, phi_etap, acc_eta, nup, 
     sig2ep, sig2etap, sig2betap, betap, betasp, rhop, mu_lp, sig2lp, w0p, wp,
     zfit);
     
     accept1 += acc_eta[0];
        
     phi_etaf[i] = phi_etap[0];
     nupf[i] = nup[0];
     
     sig2ef[i] = sig2ep[0];
     sig2etaf[i] = sig2etap[0];
     sig2betaf[i] = sig2betap[0];
     rhof[i] = rhop[0];
     for(j=0; j < r1; j++) {
         sig2lf[j+i*r1] = sig2lp[j];     
     }
     for(j=0; j < p1; j++){
         betaf[j+i*p1] = betap[j];
     }         
     for(j=0; j < m1*q1; j++){
         betaspf[j+i*m1*q1] = betasp[j];
     }         
     for(j=0; j < r1; j++) {
         mu_lf[j+i*r1] = mu_lp[j];     
     }
     for(j=0; j<m1*r1; j++){
         w0f[j+i*m1*r1] = w0p[j];
     }     
     for(j=0; j<m1*r1*T1; j++){
         wf[j+i*m1*r1*T1] = wp[j];
     }     

// for pmcc, fitted sum values    
     for(j=0; j < N1; j++){
         if(i >= brin){  
          mn_rep[j] += zfit[j];
          var_rep[j] += zfit[j]*zfit[j];
         } 
     }

       ext_sige(phi_etap, phi_eta1);
       ext_sige(nup, nu);       
       ext_sige(sig2ep, sig2e1);
       ext_sigeta(sig2etap, sig2eta1); 
       ext_rho(rhop, rho1); 
       ext_beta(p, betap, beta1);              
//       for(j=0; j < m1; j++){
//         betas[j] = betasp[j];
//       }         

     if(cov[0]==4){
     GPPsp_para_printRnu(i, its1, rep1, p1, accept1, phi_etap, nup, rhop, sig2ep, sig2etap, sig2betap, betap); 
     }                   
     else {              
     GPPsp_para_printR(i, its1, rep1, p1, accept1, phi_etap, rhop, sig2ep, sig2etap, sig2betap, betap); 
     }

    
     } // end of iteration loop
     PutRNGstate();

     free(phi_etap); free(nu); free(nup); free(sig2ep); free(sig2etap); free(sig2betap);
     free(rhop); free(betap); free(betasp); free(mu_lp); free(sig2lp); free(wp); free(w0p); 
     free(zfit); free(phi_eta1); free(sig2e1); free(sig2eta1); free(sig2beta1); free(rho1); 
     free(beta1); free(mu_l1); free(sig2l1); free(acc_eta); free(out); 

     accept_etaf[0] = accept1;

	int *iit;
    iit = (int *) malloc((size_t)((col)*sizeof(int)));
	iit[0] = its[0] - burnin[0];

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit[0];
          var_rep[j] = var_rep[j]/iit[0];
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          z_mean_sd[j] = mn_rep[j];
          z_mean_sd[N1+j] = sqrt(var_rep[j]);
     }

     double pen, go;
     pen = 0;
     go =0;

// pmcc          
     for(j=0; j < N1; j++){
         if (flag[j] == 1.0){
          mn_rep[j] = 0.0;
          var_rep[j] = 0.0;
         }
         else{
          mn_rep[j] = mn_rep[j];
          var_rep[j] = var_rep[j];
          mn_rep[j] = (mn_rep[j] - z[j])*(mn_rep[j] - z[j]);          
         }
         pen += var_rep[j]; 
         go += mn_rep[j];
     }
     penalty[0] = pen;
     gof[0] = go;       

     free(mn_rep); free(var_rep); free(iit);

     return;
}


/////////////////////////// THE END ////////////////////////////
