/**********************************************************************
 * Given nf left and right vectors such that the first ni-1 vectors are
 * bi-orthogonal, biortho_local extends the biorthogonality through the 
 * full set using Gram-Shmidt procedure. 
 * VL (IN/OUT) array of dimension ldvl*nf of the left eigenvectors stored 
 *             column by column such that VL[i*ldvl] is the first element
 *             of the left (i+1)th std::vector.
 * VR (IN/OUT) array of dimension ldvr*nf of the right eigenvectors stored 
 *             column by column such that VR[i*ldvr] is the first element
 *             of the right (i+1)th std::vector.
 *
 * biortho_global does similar things for single precision complex vectors
 * that are distributed over nodes.
 ****************************************************************************/





#include "qdp-lapack_biortho.h"

void biortho_local(Complex_Z *VL, int ldvl, Complex_Z *VR, int ldvr, int m, int ni, int nf, int ktimes){

    int ONE=1, i,j,k,iq,ielem;
    Complex_Z z_one ={+1.0e+00,+0.0e+00}, z_zero={+0.0e+00,+0.0e+00};
    Complex_Z z_mone={-1.0e+00,+0.0e+00};
    Complex_Z tmpz, tmpz1,beta;
    double tmpd,delta;


    for(i = ni-1; i <nf; i++){

      for(k=1; k <= ktimes; k++){

        
         for(j=0; j <i; j++)
            {
               wrap_zdot_small(&tmpz,&m, &VL[j*ldvl],&ONE,&VR[i*ldvr],&ONE);
               for(iq=0; iq<m; iq++)
                 {
                   z_mul_primme(&tmpz1,&tmpz,&VR[j*ldvr+iq]);
                   VR[i*ldvr+iq].r = VR[i*ldvr+iq].r - tmpz1.r;
                   VR[i*ldvr+iq].i = VR[i*ldvr+iq].i - tmpz1.i;
                 }

               
               wrap_zdot_small(&tmpz,&m, &VR[j*ldvr], &ONE, &VL[i*ldvl],&ONE);
               for(iq=0; iq<m; iq++)
                 {
                   z_mul_primme(&tmpz1,&tmpz,&VL[j*ldvl+iq]);
                   VL[i*ldvl+iq].r = VL[i*ldvl+iq].r - tmpz1.r;
                   VL[i*ldvl+iq].i = VL[i*ldvl+iq].i - tmpz1.i;
                 }

            }

         
         /* Normalize the new vectors */
         wrap_zdot_small(&tmpz,&m, &VR[i*ldvr], &ONE, &VR[i*ldvr], &ONE);
         tmpz.r = sqrt(tmpz.r);  tmpz.i=0.0;
         z_div_primme(&tmpz1, &z_one, &tmpz);  /* 1/tempz */
         BLAS_ZSCAL(&m, &tmpz1, &VR[i*ldvr], &ONE); 

         wrap_zdot_small(&tmpz,&m, &VL[i*ldvl], &ONE, &VL[i*ldvl], &ONE);
         tmpz.r=sqrt(tmpz.r); tmpz.i=0.0;
         z_div_primme(&tmpz1, &z_one, &tmpz);  /* 1/tempz */
         BLAS_ZSCAL(&m, &tmpz1, &VL[i*ldvl], &ONE);

         /* bi-orthogonalize */
         wrap_zdot_small(&tmpz,&m, &VR[i*ldvr], &ONE, &VL[i*ldvl], &ONE);
         z_div_primme(&tmpz1,&z_one,&tmpz);
         BLAS_ZSCAL(&m, &tmpz1, &VL[i*ldvl], &ONE);
         
         //tmpd   = z_abs_primme(tmpz);
         //delta  = sqrt(tmpd);
         //beta.r = tmpz.r/delta; beta.i=tmpz.i/delta;
         
         //z_div_primme(&tmpz1, &z_one, &beta);  /* 1/tempz */
         //BLAS_ZSCAL(&m, &tmpz1, &VL[i*ldvl], &ONE);

         //for(ielem=0; ielem<m; ielem++){
        //     VR[i*ldvr+ielem].r =  VR[i*ldvr+ielem].r/delta;
        //     VR[i*ldvr+ielem].i =  VR[i*ldvr+ielem].i/delta;}


      }
      
   }


 }

/*----------------------------------------------------------------------------------------------------*/
void biortho_global_C(Complex_C *VL, int ldvl, Complex_C *VR, int ldvr, int m, int ni, int nf, int ktimes, void *params){

    int ONE=1, i,j,k,iq,ielem;
    Complex_C c_one={+1.0e+00,+0.0e+00}, c_zero={+0.0e+00,+0.0e+00};
    Complex_C c_mone={-1.0e+00,+0.0e+00};
    Complex_C tmpc, tmpc1,tmpc2,beta;
    float tmpf,delta;


    for(i = ni-1; i <nf; i++){

      for(k=1; k <= ktimes; k++){

        
         for(j=0; j <i; j++)
            {
               tmpc = wrap_zsum_cdot(&m,&VL[j*ldvl],&ONE,&VR[i*ldvr],&ONE,params);
               for(iq=0; iq<m; iq++)
                 {
                   c_mul_primme(&tmpc1,&tmpc,&VR[j*ldvr+iq]);
                   VR[i*ldvr+iq].r = VR[i*ldvr+iq].r - tmpc1.r;
                   VR[i*ldvr+iq].i = VR[i*ldvr+iq].i - tmpc1.i;
                 }

               
               
               tmpc = wrap_zsum_cdot(&m,&VR[j*ldvr],&ONE,&VL[i*ldvl],&ONE,params);
               for(iq=0; iq<m; iq++)
                 {
                   c_mul_primme(&tmpc1,&tmpc,&VL[j*ldvl+iq]);
                   VL[i*ldvl+iq].r = VL[i*ldvl+iq].r - tmpc1.r;
                   VL[i*ldvl+iq].i = VL[i*ldvl+iq].i - tmpc1.i;
                 }

            }

         
         /* Normalize the new vectors */
         
         tmpc = wrap_zsum_cdot(&m,&VR[i*ldvr],&ONE,&VR[i*ldvr],&ONE,params);
         tmpc.r = sqrt(tmpc.r);  tmpc.i=0.0;
         c_div_primme(&tmpc1, &c_one, &tmpc);  
         BLAS_CSCAL(&m, &tmpc1, &VR[i*ldvr], &ONE);  

         tmpc = wrap_zsum_cdot(&m, &VL[i*ldvl],&ONE,&VL[i*ldvl],&ONE,params);
         tmpc.r=sqrt(tmpc.r); tmpc.i=0.0;
         c_div_primme(&tmpc1, &c_one, &tmpc);  
         BLAS_CSCAL(&m, &tmpc1, &VL[i*ldvl], &ONE); 

         /* bi-orthogonalize */
         tmpc = wrap_zsum_cdot(&m,&VR[i*ldvr],&ONE,&VL[i*ldvl],&ONE,params);
         c_div_primme(&tmpc1, &c_one, &tmpc);
         BLAS_CSCAL(&m, &tmpc1, &VL[i*ldvl], &ONE);

         //tmpf = c_abs_primme(tmpc);
         //delta= sqrt(tmpf);
         //beta.r = tmpc.r/delta; beta.i = tmpc.i/delta;


         //c_div_primme(&tmpc1, &c_one, &beta);  /* 1/tempc */
         //BLAS_CSCAL(&m, &tmpc1, &VL[i*ldvl], &ONE);

         //for(ielem=0; ielem<m; ielem++){
         //    VR[i*ldvr+ielem].r = VR[i*ldvr+ielem].r/delta;
         //    VR[i*ldvr+ielem].i = VR[i*ldvr+ielem].i/delta;}

      }
      
   }


 }

/*----------------------------------------------------------------------------------------------------*/


void biortho_global_Z(Complex_Z *VL, int ldvl, Complex_Z *VR, int ldvr, int m, int ni, int nf, int ktimes, void *params){

    int ONE=1, i,j,k,iq,ielem;
    Complex_Z c_one={+1.0e+00,+0.0e+00}, c_zero={+0.0e+00,+0.0e+00};
    Complex_Z c_mone={-1.0e+00,+0.0e+00};
    Complex_Z tmpc, tmpc1,tmpc2,beta;
    double tmpf,delta;


    for(i = ni-1; i <nf; i++){

      for(k=1; k <= ktimes; k++){

        
         for(j=0; j <i; j++)
            {
               tmpc = wrap_zdot(&m,&VL[j*ldvl],&ONE,&VR[i*ldvr],&ONE,params);
               for(iq=0; iq<m; iq++)
                 {
                   z_mul_primme(&tmpc1,&tmpc,&VR[j*ldvr+iq]);
                   VR[i*ldvr+iq].r = VR[i*ldvr+iq].r - tmpc1.r;
                   VR[i*ldvr+iq].i = VR[i*ldvr+iq].i - tmpc1.i;
                 }

               
               
               tmpc = wrap_zdot(&m,&VR[j*ldvr],&ONE,&VL[i*ldvl],&ONE,params);
               for(iq=0; iq<m; iq++)
                 {
                   z_mul_primme(&tmpc1,&tmpc,&VL[j*ldvl+iq]);
                   VL[i*ldvl+iq].r = VL[i*ldvl+iq].r - tmpc1.r;
                   VL[i*ldvl+iq].i = VL[i*ldvl+iq].i - tmpc1.i;
                 }

            }

         
         /* Normalize the new vectors */
         
         tmpc = wrap_zdot(&m,&VR[i*ldvr],&ONE,&VR[i*ldvr],&ONE,params);
         tmpc.r = sqrt(tmpc.r);  tmpc.i=0.0;
         z_div_primme(&tmpc1, &c_one, &tmpc);  
         BLAS_ZSCAL(&m, &tmpc1, &VR[i*ldvr], &ONE);  

         tmpc = wrap_zdot(&m, &VL[i*ldvl],&ONE,&VL[i*ldvl],&ONE,params);
         tmpc.r=sqrt(tmpc.r); tmpc.i=0.0;
         z_div_primme(&tmpc1, &c_one, &tmpc);  
         BLAS_ZSCAL(&m, &tmpc1, &VL[i*ldvl], &ONE); 

         /* bi-orthogonalize */
         tmpc = wrap_zdot(&m,&VR[i*ldvr],&ONE,&VL[i*ldvl],&ONE,params);
         z_div_primme(&tmpc1, &c_one, &tmpc);
         BLAS_ZSCAL(&m, &tmpc1, &VL[i*ldvl], &ONE);

         //tmpf = c_abs_primme(tmpc);
         //delta= sqrt(tmpf);
         //beta.r = tmpc.r/delta; beta.i = tmpc.i/delta;


         //c_div_primme(&tmpc1, &c_one, &beta);  /* 1/tempc */
         //BLAS_CSCAL(&m, &tmpc1, &VL[i*ldvl], &ONE);

         //for(ielem=0; ielem<m; ielem++){
         //    VR[i*ldvr+ielem].r = VR[i*ldvr+ielem].r/delta;
         //    VR[i*ldvr+ielem].i = VR[i*ldvr+ielem].i/delta;}

      }
      
   }


 }






