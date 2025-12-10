/*-----------------------------------------------------------------------------
  Computing and sorting eignvalues and eignvectors for a general complex matrix.
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
------------------------------------------------------------------------------
*/

#include "G_eval.h"

void ZG_eval(Complex_Z *A, int N, int LDA, Complex_Z *W, char SRT_OPT, double epsi, 
             Complex_Z *VL, int LDVL, Complex_Z *VR, int LDVR, int *info)
{
     char cV='V',cN='N';
     Complex_Z *WORK,*tmpVL,*tmpVR,*tmpA,tmpz1,tmpz2,z_one,z_zero,tmpz3;
     double    *RWORK,tmpd;
     int LWORK,i,j,k,ONE,allelems,itmp,*IPIV,info_zgesv;
     unsigned long int *indx;
     IPIV = (int *) calloc(N,sizeof(int));
     tmpA = (Complex_Z *) calloc(N*LDA,sizeof(Complex_Z));
     
  


     ONE=1;
     z_one.r=1.0; z_one.i=0.0;
     z_zero.r=0.0;  z_zero.i=0.0;



     allelems = LDA*N;
     BLAS_ZCOPY(&allelems,A,&ONE,tmpA,&ONE);   


     LWORK=3*N;

     WORK=(Complex_Z *) calloc(LWORK,sizeof(Complex_Z));
     if(WORK == NULL){
       fprintf(stderr,"Error: couldn't allocate WROK array inside ZG_eval.\n");
       exit(0);}

     RWORK=(double *) calloc(2*N,sizeof(double));
     if(RWORK == NULL){
        fprintf(stderr,"Error: couldn't allocate RWORK array inside ZG_eval.\n");
        exit(0);}


     BLAS_ZGEEV(&cN,&cV,&N,A,&LDA,W,VL,&LDVL,VR,&LDVR,WORK,&LWORK,RWORK,info);

     /* sort eigenvalues and corresponding eigenvectors if needed */
     if( (SRT_OPT == 'R') || (SRT_OPT == 'M') )
         {
            indx=(unsigned long int *) calloc(N,sizeof(unsigned long int));
            if(indx == NULL){
               fprintf(stderr,"Error: couldn't allocate indx array inside ZG_eval.\n");
               exit(0);}

            for(i=0; i<N; i++)
                 WORK[i]=W[i];   //copy the eigenvalues into WORK

            if( SRT_OPT == 'R')
                {
                   for(i=0; i<N; i++)
                      RWORK[i]=W[i].r;   //RWORK stores the real parts of the eigenvalues
                }

            if( SRT_OPT == 'M')
                {
                   for(i=0; i<N; i++)
                      RWORK[i]=z_abs_primme(W[i]);  //RWORK stores the absolute value of the eigenvalues
                }

            index_bubble_sort(N,RWORK,indx);   
           //index sort RWORK in ascending order. RWORK is unchanged but the correct order is given by indx
            
            i = 1;            //further ordering of the conjugate pairs such that the pair memeber with +imag comes first
            while ( i<N )
              {
                  if ( Zconjugate(W[indx[i-1]],W[indx[i]],epsi) )
                    {
                      if ( W[indx[i-1]].i < 0 )
                        {
                           itmp   = indx[i-1];
                           indx[i-1] = indx[i];
                           indx[i] = itmp;
                        }
                      i = i + 2;
                    }
                  else
                      i = i + 1;
              }

            //Finally return sorted eigenvalues
            for(i=0; i<N; i++)
               W[i]=WORK[indx[i]];

            //tmpVL=(Complex_Z *) calloc(LDVL*N,sizeof(Complex_Z));
            //if( tmpVL == NULL){
            //   fprintf(stderr,"Error: couldn't allocate tmpVL inside ZG_eval.\n");
            //   exit(0);}


 
            //allelems=LDVL*N;
            //BLAS_ZCOPY(&allelems,VL,&ONE,tmpVL,&ONE);

            //for(i=0; i<N; i++)
            //  BLAS_ZCOPY(&LDVL,&tmpVL[LDVL*indx[i]],&ONE,&VL[LDVL*i],&ONE);

            //free(tmpVL);

            tmpVR=(Complex_Z *) calloc(LDVR*N,sizeof(Complex_Z));
            if(tmpVR==NULL){
               fprintf(stderr,"Error: couldn't allocate tmpVR inside ZG_eval.\n");
               exit(0);}            



            allelems=LDVR*N;
            BLAS_ZCOPY(&allelems,VR,&ONE,tmpVR,&ONE);

            for(i=0; i<N; i++)
              BLAS_ZCOPY(&N,&tmpVR[LDVR*indx[i]],&ONE,&VR[LDVR*i],&ONE);


            //Check the computation of the right eignvectors
            /*
            for(i=0; i<N; i++)
               {
                   BLAS_ZCOPY(&N,&VR[LDVR*i],&ONE,&WORK[N],&ONE);
                   BLAS_ZGEMV(&cN,&N,&N,&z_one,tmpA,&LDA,&VR[LDVR*i],&ONE,&z_zero,WORK,&ONE);  
                   BLAS_ZSCAL(&N,&W[i],&WORK[N],&ONE);
                   for(j=0; j<N; j++)
                     {  WORK[j].r = WORK[N+j].r - WORK[j].r;
                        WORK[j].i = WORK[N+j].i - WORK[j].i;
                     }
                   wrap_zdot_small(&tmpz1,&N,WORK,&ONE,WORK,&ONE);
                   printf("Eval# %d, value %g %g , resnorm= %g\n",i,W[i].r,W[i].i,sqrt(tmpz1.r));
              }
             */








            for(i=0; i<N; i++)
               for(j=0; j<N; j++)
                 {
                     tmpVR[i+j*LDVR].r  = VR[j+i*LDVR].r;
                     tmpVR[i+j*LDVR].i  = -VR[j+i*LDVR].i;
                 }


            //define a unit matrix in VL
            for(i=0; i<LDVL*N; i++)
               {VL[i].r=0.0; VL[i].i=0.0;}


            for(i=0; i<N; i++)
                VL[i+i*LDVL].r=1.0;

            BLAS_ZGESV(&N,&N,tmpVR,&LDVR,IPIV,VL,&LDVL,&info_zgesv);
            if(info_zgesv != 0)
               {printf("error solving for left eigenvectors in zg_eval.\n");
                exit(0);}

            //check biorthogonality of the left and right vectors
            //printf("Biorthogonality of VL and VR in ZG_eval:\n");
            //for(i=0; i<N; i++)
            //  for(j=0; j<=i; j++)
            //     {
            //        wrap_zdot_small(&tmpz1,&N,&VL[i*LDVL],&ONE,&VR[j*LDVR],&ONE);
            //        printf("%d %d %g %g\n",i,j,tmpz1.r,tmpz1.i);
            //     }

         }  //end if( (SRT_OPT == 'R') || (SRT_OPT == 'M') )
     else
       {
           fprintf(stderr,"Invalid choice of SRT_OPT. EIgenvalues are not sorted\n");
           exit(0);
        }


     /* Normalize eigenvectors such that VL(:,k)^H*VR(:,k) = 1
        it is assumed that the matrix is non-defective, i.e. 
        VL(:,k)^H * VR(:,j) =0.0. If this is not the case then
        further consideration is needed. */

     /*
     for(i=0; i<N; i++)
        {
           wrap_zdot_small(&tmpz1,&N,&VL[LDVL*i],&ONE,&VL[LDVL*i],&ONE);
           tmpz1.r=sqrt(fabs((double) tmpz1.r));
           tmpz1.i=0.0;
           z_div_primme(&tmpz2,&z_one,&tmpz1);
           BLAS_ZSCAL(&N,&tmpz2,&VL[LDVL*i],&ONE);

           wrap_zdot_small(&tmpz1,&N,&VR[LDVR*i],&ONE,&VR[LDVR*i],&ONE);
           tmpz1.r=sqrt(fabs((double) tmpz1.r));
           tmpz1.i=0.0;
           z_div_primme(&tmpz2,&z_one,&tmpz1);
           BLAS_ZSCAL(&N,&tmpz2,&VR[LDVR*i],&ONE);

           wrap_zdot_small(&tmpz1,&N,&VR[LDVR*i],&ONE,&VL[LDVL*i],&ONE);

           tmpd = z_abs_primme(tmpz1);
           tmpd = sqrt(tmpd);
           tmpz3.r = tmpz1.r/tmpd; tmpz3.i = tmpz1.i/tmpd;
           z_div_primme(&tmpz2,&z_one,&tmpz3);
           BLAS_ZSCAL(&N,&tmpz2,&VL[LDVL*i],&ONE);

           for(k=0; k<N; k++){
              VR[LDVR*i+k].r = VR[LDVR*i+k].r/tmpd;
              VR[LDVR*i+k].i = VR[LDVR*i+k].i/tmpd;}

        }
      */
}



/***************************************************************************************************/
void ZG_eval_original(Complex_Z *A, int N, int LDA, Complex_Z *W, char SRT_OPT, double epsi, 
             Complex_Z *VL, int LDVL, Complex_Z *VR, int LDVR, int *info)
{
     char cV='V';
     Complex_Z *WORK,*tmpVL,*tmpVR,tmpz1,tmpz2,z_one,tmpz3;
     double    *RWORK,tmpd;
     int LWORK,i,j,k,ONE,allelems,itmp;
     unsigned long int *indx;
     
     ONE=1;
     z_one.r=1.0; z_one.i=0.0;

     LWORK=3*N;

     WORK=(Complex_Z *) calloc(LWORK,sizeof(Complex_Z));
     if(WORK == NULL){
       fprintf(stderr,"Error: couldn't allocate WROK array inside ZG_eval.\n");
       exit(0);}

     RWORK=(double *) calloc(2*N,sizeof(double));
     if(RWORK == NULL){
        fprintf(stderr,"Error: couldn't allocate RWORK array inside ZG_eval.\n");
        exit(0);}


     BLAS_ZGEEV(&cV,&cV,&N,A,&LDA,W,VL,&LDVL,VR,&LDVR,WORK,&LWORK,RWORK,info);

     /* sort eigenvalues and corresponding eigenvectors if needed */
     if( (SRT_OPT == 'R') || (SRT_OPT == 'M') )
         {
            indx=(unsigned long int *) calloc(N,sizeof(unsigned long int));
            if(indx == NULL){
               fprintf(stderr,"Error: couldn't allocate indx array inside ZG_eval.\n");
               exit(0);}

            for(i=0; i<N; i++)
                 WORK[i]=W[i];   //copy the eigenvalues into WORK

            if( SRT_OPT == 'R')
                {
                   for(i=0; i<N; i++)
                      RWORK[i]=W[i].r;   //RWORK stores the real parts of the eigenvalues
                }

            if( SRT_OPT == 'M')
                {
                   for(i=0; i<N; i++)
                      RWORK[i]=z_abs_primme(W[i]);  //RWORK stores the absolute value of the eigenvalues
                }

            index_bubble_sort(N,RWORK,indx);   
           //index sort RWORK in ascending order. RWORK is unchanged but the correct order is given by indx
            
            i = 1;            //further ordering of the conjugate pairs such that the pair memeber with +imag comes first
            while ( i<N )
              {
                  if ( Zconjugate(W[indx[i-1]],W[indx[i]],epsi) )
                    {
                      if ( W[indx[i-1]].i < 0 )
                        {
                           itmp   = indx[i-1];
                           indx[i-1] = indx[i];
                           indx[i] = itmp;
                        }
                      i = i + 2;
                    }
                  else
                      i = i + 1;
              }

            //Finally return sorted eigenvalues
            for(i=0; i<N; i++)
               W[i]=WORK[indx[i]];

            tmpVL=(Complex_Z *) calloc(LDVL*N,sizeof(Complex_Z));
            if( tmpVL == NULL){
               fprintf(stderr,"Error: couldn't allocate tmpVL inside ZG_eval.\n");
               exit(0);}


 
            allelems=LDVL*N;
            BLAS_ZCOPY(&allelems,VL,&ONE,tmpVL,&ONE);

            for(i=0; i<N; i++)
              BLAS_ZCOPY(&LDVL,&tmpVL[LDVL*indx[i]],&ONE,&VL[LDVL*i],&ONE);

            free(tmpVL);

            tmpVR=(Complex_Z *) calloc(LDVR*N,sizeof(Complex_Z));
            if(tmpVR==NULL){
               fprintf(stderr,"Error: couldn't allocate tmpVR inside ZG_eval.\n");
               exit(0);}            



            allelems=LDVR*N;
            BLAS_ZCOPY(&allelems,VR,&ONE,tmpVR,&ONE);

            for(i=0; i<N; i++)
              BLAS_ZCOPY(&LDVR,&tmpVR[LDVR*indx[i]],&ONE,&VR[LDVR*i],&ONE);

            free(tmpVR);

         }  //end if( (SRT_OPT == 'R') || (SRT_OPT == 'M') )
     else
       {
           fprintf(stderr,"Invalid choice of SRT_OPT. EIgenvalues are not sorted\n");
           exit(0);
        }


     /* Normalize eigenvectors such that VL(:,k)^H*VR(:,k) = 1
        it is assumed that the matrix is non-defective, i.e. 
        VL(:,k)^H * VR(:,j) =0.0. If this is not the case then
        further consideration is needed. */


     for(i=0; i<N; i++)
        {
           wrap_zdot_small(&tmpz1,&N,&VL[LDVL*i],&ONE,&VL[LDVL*i],&ONE);
           tmpz1.r=sqrt(fabs((double) tmpz1.r));
           tmpz1.i=0.0;
           z_div_primme(&tmpz2,&z_one,&tmpz1);
           BLAS_ZSCAL(&N,&tmpz2,&VL[LDVL*i],&ONE);

           wrap_zdot_small(&tmpz1,&N,&VR[LDVR*i],&ONE,&VR[LDVR*i],&ONE);
           tmpz1.r=sqrt(fabs((double) tmpz1.r));
           tmpz1.i=0.0;
           z_div_primme(&tmpz2,&z_one,&tmpz1);
           BLAS_ZSCAL(&N,&tmpz2,&VR[LDVR*i],&ONE);

           wrap_zdot_small(&tmpz1,&N,&VR[LDVR*i],&ONE,&VL[LDVL*i],&ONE);

           tmpd = z_abs_primme(tmpz1);
           tmpd = sqrt(tmpd);
           tmpz3.r = tmpz1.r/tmpd; tmpz3.i = tmpz1.i/tmpd;
           z_div_primme(&tmpz2,&z_one,&tmpz3);
           BLAS_ZSCAL(&N,&tmpz2,&VL[LDVL*i],&ONE);

           for(k=0; k<N; k++){
              VR[LDVR*i+k].r = VR[LDVR*i+k].r/tmpd;
              VR[LDVR*i+k].i = VR[LDVR*i+k].i/tmpd;}

        }

}

/*-----------------------------------------------------------------------------------*/


void CG_eval(Complex_C *A, int N, int LDA, Complex_C *W, char SRT_OPT, float epsi, 
             Complex_C *VL, int LDVL, Complex_C *VR, int LDVR, int *info)
{
     char cV='V';
     Complex_C *WORK,*tmpVL,*tmpVR,tmpz1,tmpz2,tmpz3,z_one;
     float    *RWORK,tmpd;
     int LWORK,i,j,k,ONE,allelems,itmp;
     unsigned long int *indx;
     
     ONE=1;
     z_one.r=1.0; z_one.i=0.0;

     LWORK=3*N;

     WORK=(Complex_C *) calloc(LWORK,sizeof(Complex_C));
     if(WORK == NULL){
       fprintf(stderr,"Error: couldn't allocate WROK array inside CG_eval.\n");
       exit(0);}

     RWORK=(float *) calloc(2*N,sizeof(float));
     if(RWORK == NULL){
       fprintf(stderr,"Error: couldn't allocate RWROK array inside CG_eval.\n");
       exit(0);}



     BLAS_CGEEV(&cV,&cV,&N,A,&LDA,W,VL,&LDVL,VR,&LDVR,WORK,&LWORK,RWORK,info);

     /* sort eigenvalues and corresponding eigenvectors if needed */
     if( (SRT_OPT == 'R') || (SRT_OPT == 'M') )
         {
            indx=(unsigned long int *) calloc(N,sizeof(unsigned long int));
            if( indx == NULL){
               fprintf(stderr,"Error: couldn't allocate indx array inside CG_eval.\n");
               exit(0);}

            for(i=0; i<N; i++)
                      WORK[i]=W[i];

            if( SRT_OPT == 'R')
                {
                   for(i=0; i<N; i++)
                      RWORK[i]=W[i].r;
                }

            if( SRT_OPT == 'M')
                {
                   for(i=0; i<N; i++)
                      RWORK[i]=c_abs_primme(W[i]);
                }

            index_bubble_sort_single(N,RWORK,indx);
            
            i = 1;   
            while ( i<N )
              {
                  if ( Cconjugate(W[indx[i-1]],W[indx[i]],epsi) )
                    {
                      if ( W[indx[i-1]].i < 0 )
                        {
                           itmp   = indx[i-1];
                           indx[i-1] = indx[i];
                           indx[i] = itmp;
                        }
                      i = i + 2;
                    }
                  else
                      i = i + 1;
              }

            for(i=0; i<N; i++)
               W[i]=WORK[indx[i]];

            tmpVL=(Complex_C *) calloc(LDVL*N,sizeof(Complex_C));
            if(tmpVL == NULL){
               fprintf(stderr,"Error: couldn't allocate tmpVL array inside CG_eval.\n");
               exit(0);}

 
            allelems=LDVL*N;
            BLAS_CCOPY(&allelems,VL,&ONE,tmpVL,&ONE);
            for(i=0; i<N; i++)
              BLAS_CCOPY(&LDVL,&tmpVL[LDVL*indx[i]],&ONE,&VL[LDVL*i],&ONE);
            free(tmpVL);

            tmpVR=(Complex_C *) calloc(LDVR*N,sizeof(Complex_C));
            if( tmpVR == NULL){
               fprintf(stderr,"Error: couldn't allocate tmpVR array inside CG_eval.\n");
               exit(0);}



            allelems=LDVR*N;
            BLAS_CCOPY(&allelems,VR,&ONE,tmpVR,&ONE);
            for(i=0; i<N; i++)
              BLAS_CCOPY(&LDVR,&tmpVR[LDVR*indx[i]],&ONE,&VR[LDVR*i],&ONE);
            free(tmpVR);
         }
     else
         {
             printf("Invalid choice of SRT_OPT. EIgenvalues are not sorted\n");
             exit(0);
         }

     /* Normalize eigenvectors such that VL(:,k)^H*VR(:,k) = 1
        it is assumed that the matrix is non-defective, i.e. 
        VL(:,k)^H * VR(:,j) =0.0. If this is not the case then
        further consideration is needed. */
     for(i=0; i<N; i++)
        {
           wrap_cdot_small(&tmpz1,&N,&VL[LDVL*i],&ONE,&VL[LDVL*i],&ONE);
           tmpz1.r=sqrt(tmpz1.r);
           tmpz1.i=0.0;
           c_div_primme(&tmpz2,&z_one,&tmpz1);
           BLAS_CSCAL(&N,&tmpz2,&VL[LDVL*i],&ONE);

           wrap_cdot_small(&tmpz1,&N,&VR[LDVR*i],&ONE,&VR[LDVR*i],&ONE);
           tmpz1.r=sqrt(tmpz1.r);
           tmpz1.i=0.0;
           c_div_primme(&tmpz2,&z_one,&tmpz1);
           BLAS_CSCAL(&N,&tmpz2,&VR[LDVR*i],&ONE);

           wrap_cdot_small(&tmpz1,&N,&VR[LDVR*i],&ONE,&VL[LDVL*i],&ONE);
           tmpd = c_abs_primme(tmpz1);
           tmpd = sqrt(tmpd);
           tmpz3.r = tmpz1.r/tmpd; tmpz3.i = tmpz1.i/tmpd;
           c_div_primme(&tmpz2,&z_one,&tmpz3);
           BLAS_CSCAL(&N,&tmpz2,&VL[LDVL*i],&ONE);

           for(k=0; k<N; k++){
              VR[LDVR*i+k].r = VR[LDVR*i+k].r/tmpd;
              VR[LDVR*i+k].i = VR[LDVR*i+k].i/tmpd;}
        }

}

/*-----------------------------------------------------------------------------------*/

