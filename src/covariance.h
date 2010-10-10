
#include "header.h"

/*****************************************************************************/
/*****************************************************************************/

/****************************** From "covariance.c" file ***************************/
/*****************************************************************************/

void covFormat(int *cov, int *n, double *phi, double *d, double *sig2eta, 
     double *S, double *det, double *Sinv, double *Qeta);
void covFormat1(int *cov, int *n, double *phi, double *d, double *S);
void covFormat2(int *cov, int *n, double *phi, double *d, double *det, 
     double *Sinv);
void covFormat3(int *cov, int *n, double *phi, double *d, double *sig_eta, 
     double *det, double *Qeta);
     
void covExpo(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta);
void covGaus(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta);
void covSphe(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta);
void covMatern32(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta);
void covMatern(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta);
     
void covExpo1(int *n, double *phi, double *d, double *S);
void covGaus1(int *n, double *phi, double *d, double *S);
void covSphe1(int *n, double *phi, double *d, double *S);
void covMatern321(int *n, double *phi, double *d, double *S);
void covMatern1(int *n, double *phi, double *d, double *S);
     
void covExpo2(int *n, double *phi, double *d, double *det, double *Sinv);
void covGaus2(int *n, double *phi, double *d, double *det, double *Sinv);
void covSphe2(int *n, double *phi, double *d, double *det, double *Sinv);
void covMatern322(int *n, double *phi, double *d, double *det, double *Sinv);
void covMatern2(int *n, double *phi, double *d, double *det, double *Sinv);

double bessk1(double x);
double bessi1(double x);
          
void Q_eta(double *phi, double *sig_eta, int *n, double *d, double *Qeta);
void S_inv(double *phi, double *d, int *n, int *constant, double *Sinv);
void S_inv_det(double *phi, double *d, int *n, int *constant, double *det, 
     double *Sinv);


     
/*****************************************************************************/
/*****************************************************************************/
