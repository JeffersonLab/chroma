/*------------------------------------------------------------------------
  BICGSTAB algorithm for solving A*x=b with A a sparse non-symmetric matrix 
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
--------------------------------------------------------------------------*/

#include "bicgstab.h"

void MyBICGSTAB_C(int n, int lde, Complex_C *x, Complex_C *b, float tol, int maxiter, int *iter, float *res,
                   void (*matvec)(void *, void *, void *), void *params, float AnormEst, int ConvTestOpt,
                   Complex_C *work, int *flag)
{

     Complex_C *r,*rs,*p,*ap,*s,*as;
     Complex_C theta,theta_prev,alpha,beta,omega;
     Complex_C tmpc,tmpc1,tmpc2,tmpc3;
     Complex_C tzero={0.0e+00,0.0e+00};
     Complex_C tpone={1.0e+00,0.0e+00};
     int i,j,k,ONE=1;
     float stoptol_used,stoptol_cur,xnorm,bnorm;
     float MACHEPS = 1e-7;

     
     r=work; 
     rs=work+lde;
     p=work +2*lde;
     ap=work+3*lde;
     s=work +4*lde;
     as=work+5*lde;

     
     tmpc = wrap_zsum_cdot(&n,b,&ONE,b,&ONE,params);
     bnorm = sqrt(tmpc.r);
     
     

     if(bnorm == 0)
       {
           for(i=0; i<n; i++)
               x[i]=tzero;
           (*iter)=0;
           (*flag)=1;
           return ;
       }


     tmpc = wrap_zsum_cdot(&n,x,&ONE,x,&ONE,params);
     xnorm = sqrt(tmpc.r);
     if(xnorm > 0.0)
       {
           matvec(x,r,params);
           for(i=0; i<n; i++)
              {
                  r[i].r = b[i].r - r[i].r;
                  r[i].i = b[i].i - r[i].i;
              }
        }
      else
        BLAS_CCOPY(&n,b,&ONE,r,&ONE);


      tmpc = wrap_zsum_cdot(&n,r,&ONE,r,&ONE,params);
      res[0]=sqrt(tmpc.r);

      for(i=0; i<n; i++)
          rs[i]=r[i]; 
         

      j=0;
      while(j<maxiter-1)
        {
            j = j+1;
            theta = wrap_zsum_cdot(&n,rs,&ONE,r,&ONE,params);
            if(c_abs_primme(theta)==0)
               {
                  (*iter)=j;
                  (*flag)=2;
                  return;
               }
            if(j == 1)
               BLAS_CCOPY(&n,r,&ONE,p,&ONE);
            else
               {
                  c_div_primme(&tmpc1,&theta,&theta_prev);
                  c_div_primme(&tmpc2,&alpha,&omega);
                  c_mul_primme(&beta,&tmpc1,&tmpc2);
                  tmpc.r=-omega.r; tmpc.i=-omega.i;
                  BLAS_CAXPY(&n,&tmpc,ap,&ONE,p,&ONE);      
                  BLAS_CSCAL(&n,&beta,p,&ONE);
                  BLAS_CAXPY(&n,&tpone,r,&ONE,p,&ONE);
               }


             matvec(p,ap,params);

             tmpc = wrap_zsum_cdot(&n,rs,&ONE,ap,&ONE,params);
             c_div_primme(&alpha,&theta,&tmpc);
             BLAS_CCOPY(&n,r,&ONE,s,&ONE);
             tmpc.r=-alpha.r; tmpc.i=-alpha.i;
             BLAS_CAXPY(&n,&tmpc,ap,&ONE,s,&ONE);

             BLAS_CAXPY(&n,&alpha,p,&ONE,x,&ONE);
             
             matvec(s,as,params);
             tmpc = wrap_zsum_cdot(&n,as,&ONE,s,&ONE,params);
             tmpc1 = wrap_zsum_cdot(&n,as,&ONE,as,&ONE,params);
             omega.r=tmpc.r/tmpc1.r; omega.i=tmpc.i/tmpc1.r; 
             BLAS_CAXPY(&n,&omega,s,&ONE,x,&ONE);

             BLAS_CCOPY(&n,s,&ONE,r,&ONE);
             tmpc.r=-omega.r; tmpc.i=-omega.i;
             BLAS_CAXPY(&n,&tmpc,as,&ONE,r,&ONE);
             tmpc = wrap_zsum_cdot(&n,r,&ONE,r,&ONE,params);
             res[j]=sqrt(tmpc.r);

             
             stoptol_used = tol*bnorm;
             if( ConvTestOpt == 2)
               {
                   tmpc = wrap_zsum_cdot(&n,x,&ONE,x,&ONE,params);
                   xnorm= sqrt(tmpc.r);
                   stoptol_cur = MACHEPS*(AnormEst*xnorm+bnorm);
                   if( stoptol_used < stoptol_cur )
                      stoptol_used = stoptol_cur;
               }
             
             if(res[j]< stoptol_used)
              { 
                (*iter) = j;   
                (*flag) = 0;
                printf("stopped when norm(r) <=  %g \n",stoptol_used);
                return;
              }

            if(c_abs_primme(omega)==0)
               {
                  (*iter)=j;
                  (*flag)=2;
                  return;
               }
            theta_prev.r= theta.r;
            theta_prev.i= theta.i;

           
        }

      //if no convergence after maxiter
      (*iter)=j;
      (*flag)=3;
      return;
}


/*------------------------------------------------------------------------------------------------*/

void MyBICGSTAB_Z(int n, int lde, Complex_Z *x, Complex_Z *b, double tol, int maxiter, int *iter, double *res,
                   void (*matvec)(void *, void *, void *), void *params, double AnormEst, int ConvTestOpt,
                   Complex_Z *work, int *flag)
{

     Complex_Z *r,*rs,*p,*ap,*s,*as;
     Complex_Z theta,theta_prev,alpha,beta,omega;
     Complex_Z tmpc,tmpc1,tmpc2,tmpc3;
     Complex_Z tzero={0.0e+00,0.0e+00};
     Complex_Z tpone={1.0e+00,0.0e+00};
     int i,j,k,ONE=1;
     double stoptol_used,stoptol_cur,xnorm,bnorm;
     double MACHEPS = 1e-16;

     
     r=work; 
     rs=work+lde;
     p=work +2*lde;
     ap=work+3*lde;
     s=work +4*lde;
     as=work+5*lde;


     tmpc = wrap_zdot(&n,b,&ONE,b,&ONE,params);
     bnorm = sqrt(tmpc.r);
     
     if(bnorm == 0)
       {
           for(i=0; i<n; i++)
               x[i]=tzero;
           (*iter)=0;
           (*flag)=1;
           return ;
       }


     tmpc = wrap_zdot(&n,x,&ONE,x,&ONE,params);
     xnorm = sqrt(tmpc.r);
     if(xnorm > 0.0)
       {
           matvec(x,r,params);
           for(i=0; i<n; i++)
              {
                  r[i].r = b[i].r - r[i].r;
                  r[i].i = b[i].i - r[i].i;
              }
        }
      else
        BLAS_ZCOPY(&n,b,&ONE,r,&ONE);


      tmpc = wrap_zdot(&n,r,&ONE,r,&ONE,params);
      res[0]=sqrt(tmpc.r);

      for(i=0; i<n; i++)
          rs[i]=r[i]; 
         

      j=0;
      while(j<maxiter-1)
        {
            j = j+1;
            theta = wrap_zdot(&n,rs,&ONE,r,&ONE,params);
            if(z_abs_primme(theta)==0)
               {
                  (*iter)=j;
                  (*flag)=2;
                  return;
               }
            if(j == 1)
               BLAS_ZCOPY(&n,r,&ONE,p,&ONE);
            else
               {
                  z_div_primme(&tmpc1,&theta,&theta_prev);
                  z_div_primme(&tmpc2,&alpha,&omega);
                  z_mul_primme(&beta,&tmpc1,&tmpc2);
                  tmpc.r=-omega.r; tmpc.i=-omega.i;
                  BLAS_ZAXPY(&n,&tmpc,ap,&ONE,p,&ONE);      
                  BLAS_ZSCAL(&n,&beta,p,&ONE);
                  BLAS_ZAXPY(&n,&tpone,r,&ONE,p,&ONE);
               }


             matvec(p,ap,params);

             tmpc = wrap_zdot(&n,rs,&ONE,ap,&ONE,params);
             z_div_primme(&alpha,&theta,&tmpc);
             BLAS_ZCOPY(&n,r,&ONE,s,&ONE);
             tmpc.r=-alpha.r; tmpc.i=-alpha.i;
             BLAS_ZAXPY(&n,&tmpc,ap,&ONE,s,&ONE);

             BLAS_ZAXPY(&n,&alpha,p,&ONE,x,&ONE);
             
             matvec(s,as,params);
             tmpc = wrap_zdot(&n,as,&ONE,s,&ONE,params);
             tmpc1 = wrap_zdot(&n,as,&ONE,as,&ONE,params);
             omega.r=tmpc.r/tmpc1.r; omega.i=tmpc.i/tmpc1.r; 
             BLAS_ZAXPY(&n,&omega,s,&ONE,x,&ONE);
             BLAS_ZCOPY(&n,s,&ONE,r,&ONE);
             tmpc.r=-omega.r; tmpc.i=-omega.i;
             BLAS_ZAXPY(&n,&tmpc,as,&ONE,r,&ONE);
             tmpc = wrap_zdot(&n,r,&ONE,r,&ONE,params);
             res[j]=sqrt(tmpc.r);

             
             stoptol_used = tol*bnorm;
             if( ConvTestOpt == 2)
               {
                   tmpc = wrap_zdot(&n,x,&ONE,x,&ONE,params);
                   xnorm= sqrt(tmpc.r);
                   stoptol_cur = MACHEPS*(AnormEst*xnorm+bnorm);
                   if( stoptol_used < stoptol_cur )
                      stoptol_used = stoptol_cur;
               }
             
             if(res[j]< stoptol_used)
              { 
                (*iter) = j;   
                (*flag) = 0;
                printf("stopped when norm(r) <=  %g \n",stoptol_used);
                return;
              }

            if(z_abs_primme(omega)==0)
               {
                  (*iter)=j;
                  (*flag)=2;
                  return;
               }
            theta_prev.r= theta.r;
            theta_prev.i= theta.i;

           
        }

      //if no convergence after maxiter
      (*iter)=j;
      (*flag)=3;
      return;
}


/*------------------------------------------------------------------------------------------------*/

