/*********************************************************************************
  Incremental eigbicg for solving systems with multiple right-hand sides
  Authors: Abdou Abdel-Rehim, Andreas Stathopolous, and Kostas Orignos
           andreas@cs.wm.edu
  Last updated : 10/15/2009
**********************************************************************************/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "qdp-lapack_IncrEigbicg.h"
#include "qdp-lapack_numerical.h"
#include "qdp-lapack_numerical_private.h"
#include "G_eval.h"


/*****************SINGLE PRECISION VERSION*************************************/

/******************************************************************************/
void IncrEigbicg_C( int n, int lde,int nrhs, Complex_C *X, Complex_C *B, int *ncurEvals,         	
                    int ldh, Complex_C *evecsl, Complex_C *evecsr, Complex_C *evals, 		
	            Complex_C *H,void (*matvec) (void *, void *, void *),  
	            void (*mathvec)(void *, void *, void *),void *params, float *AnormEst, 
	            Complex_C *work, Complex_C *VL,int ldvl,Complex_C *VR, int ldvr,        
	            Complex_C *ework, int esize,float tol,float *restartTol, 	
	            int maxit, char SRT_OPT, float epsi, int ConvTestOpt, int plvl,int nev,
                    int v_max,FILE *outputFile)      
{

  /* Timing vars */
  double wt1,wt2,wE,wI;

  /* Pointers */
  Complex_C tempc, tempc1, tempc2, *tmpH, *x, *resid, *b;
  float     *rnorms, *reshist, normb, curTol,resNorm,leftTol;
  int i,j,k, *IPIV, ONE = 1;
  int fs, cs, ds, tmpsize, is, phase, allelems ;
  int numIts, flag, flag2,nAdded, nev_used, iters_used, info;
  Complex_C tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
  char cR = 'R'; char cL = 'L'; char cN ='N'; 
  char cV = 'V'; char cU = 'U'; char cC ='C';
  float *xrnorms,*xlnorms,*ernorms;
  Complex_C *angles;
  

  /* ------------------------------------------------------------------- */
  /* Work allocations */
  /* ------------------------------------------------------------------- */
  fs = sizeof(float); 
  cs = sizeof(Complex_C); 
  ds = sizeof(double);
  is = sizeof(int);


  if( (IPIV = (int *) calloc(ldh,is)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate IPIV\n");
      exit(1);}


  if( (x = (Complex_C *) calloc(lde,cs)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate x\n");
      exit(1);}


  if( (resid = (Complex_C *) calloc(lde,cs)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate resid\n");
      exit(1);}

  
  if( (tmpH = (Complex_C *) calloc(ldh*ldh,cs)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate tmpH\n");
      exit(1);}


  if ((work = (Complex_C *) calloc(6*lde, cs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate work\n");
     exit(1);}
     
  if ((ework = (Complex_C *) calloc(esize, cs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate ework\n");
     exit(1);}
  
  if ((VL = (Complex_C *) calloc(v_max*ldvl, cs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate VL\n");
     exit(1);}

  if ((VR = (Complex_C *) calloc(v_max*ldvr, cs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate VR\n");
     exit(1);}
     
  if ( (rnorms = (float *) calloc(ldh, fs))  == NULL )
     {fprintf(stderr, "ERROR IncrEigbicg could not allocate rnorms\n");
      exit(1);}

  if ( (reshist = (float *) calloc(maxit, fs)) == NULL)
     {fprintf(stderr, "ERROR IncrEigbicg could not allocate reshist\n");
      exit(1);}



  if( (xlnorms = (float *) calloc(ldh,fs)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate xlnorms\n");
     exit(1);}

  if( (xrnorms = (float *) calloc(ldh,fs)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate xrnorms\n");
     exit(1);}

  if( (ernorms = (float *) calloc(ldh,fs)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate ernorms\n");
     exit(1);}

  if( (angles = (Complex_C *) calloc(ldh,cs)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate angles\n");
     exit(1);}

  /* ------------------------------------------------------------------- */
  /* end Work allocations */
  /* ------------------------------------------------------------------- */

  /* ---------------------------------------------------------------------------------  */
  /* Solving one by one the nrhs systems with incremental init-eigbicg or init-bicgstab */
  /* ---------------------------------------------------------------------------------- */


  for (j=0; j<nrhs; j++) {

      b = &B[j*lde];
      tempc = wrap_zsum_cdot(&n,b,&ONE,b,&ONE,params);
      normb = sqrt(tempc.r);

      numIts = 0;
      //choose eigbicg or bicgstab
      if(ldh-(*ncurEvals) >= nev )
        phase=1;
      else
        phase=2;

      if (plvl) fprintf(outputFile, "\n\nSystem %d\n", j);

      wE = 0.0; wI = 0.0;     /* Start accumulator timers */

      if ( (*ncurEvals > 0) && (phase==1)) {
         /* --------------------------------------------------------- */
         /* Perform init-BICG with evecsl and evecsr vectors          */
         /* xinit = xinit + evecsr*inv(H)*evecl'*(b-Ax0) 		      */
         /* --------------------------------------------------------- */
         wt1 = primme_get_wtime(); 
   
        /* copy H into tmpH otherwise it will be changed by CGEV */
        allelems = ldh*(*ncurEvals);
        BLAS_CCOPY(&allelems,H,&ONE,tmpH,&ONE);
        // LU factorization of tmpH
        BLAS_CGETRF(ncurEvals,ncurEvals,tmpH,&ldh,IPIV,&info);


        for(i=0; i<n; i++)
            x[i]=tzero;

        matvec(&X[j*lde],resid,params);

        for(i=0; i<n; i++){
         resid[i].r = b[i].r - resid[i].r;
         resid[i].i = b[i].i - resid[i].i;}


        init_BICG_C(evecsl,ldvl,evecsr,ldvr,(*ncurEvals), 
                  x,resid,lde,n,tmpH,ldh,IPIV,work,matvec,params,&info);

        if(phase==1){
           for(i=0; i<n; i++){
              X[j*lde+i].r = X[j*lde+i].r + x[i].r;
              X[j*lde+i].i = X[j*lde+i].i + x[i].i;}}
        

        wt2 = primme_get_wtime();
	wI = wI + wt2-wt1;
       }
       /* end of init-BICG with evecsl and evecsr vectors              */
       /* ------------------------------------------------------------ */

       if(phase == 1){	
         /* ------------------------------------------------------------ */
         /* Solve Ax = b with x initial guess using eigbicg and compute
            new nev eigenvectors                                         */
         /* ------------------------------------------------------------ */
         wt1 = primme_get_wtime(); 
        
         Ceigbicg(n, lde, &X[j*lde], b, &normb, tol, maxit, SRT_OPT,epsi,ConvTestOpt, &numIts, reshist, 
                 &flag, plvl, work, matvec, mathvec, params, AnormEst, nev, 
                 &evals[(*ncurEvals)], &rnorms[(*ncurEvals)],
  	         v_max, VR,ldvr,VL,ldvl,esize,ework);

         wt2 = primme_get_wtime();
	 wE = wE + wt2-wt1;

         /* ---------- */
         /* Reporting  */
         /* ---------- */
         tempc1 = wrap_zsum_cdot(&n,&X[j*lde],&ONE,&X[j*lde],&ONE,params);
         if (plvl) {
            fprintf(outputFile, "For this rhs:\n");
            fprintf(outputFile, "Norm(solution) %-16.12E, AnormEst %-16.12E\n",sqrt(tempc1.r),(*AnormEst));
            fprintf(outputFile, "Total initBICG Wallclock : %-f\n", wI);
            fprintf(outputFile, "Total eigbicg Wallclock : %-f\n", wE);
            fprintf(outputFile, "Iterations: %-d\n", numIts); 
            fprintf(outputFile, "Actual Resid of LinSys  : %e\n", reshist[numIts-1]);
	    if (plvl > 1) 
               for (i=0; i < nev; i++) 
                   fprintf(outputFile, "Eval[%d]: %-22.15E    %-22.15E           rnorm: %-22.15E\n", 
	                       i+1, evals[*ncurEvals+i].r,  evals[*ncurEvals+i].i, rnorms[*ncurEvals+i]); 
            if (plvl >1)
              {
                 fprintf(outputFile,"Residual norm\n");
               for( i=0; i < numIts; i++)
                   fprintf(outputFile,"%-d  %-22.15E\n",i,reshist[i]);
              }  

            if (flag != 0) {
               fprintf(outputFile, "Error: eigbicg returned with nonzero exit status\n");
            return;}
         }

      /* ------------------------------------------------------------------- */
      /* Update the evecsl, evecsr,  and evecsl'*A*evecsr                               */
      /* ------------------------------------------------------------------- */
      wt1 = primme_get_wtime(); 


      /* Append new Ritz pairs to evecs */
      for (i=0; i<nev; i++){
	 BLAS_CCOPY(&n, &VL[i*ldvl], &ONE, &evecsl[((*ncurEvals)+i)*ldvl], &ONE);
	 BLAS_CCOPY(&n, &VR[i*ldvr], &ONE, &evecsr[((*ncurEvals)+i)*ldvr], &ONE);}
  
      /* Bi-Orthogonalize the new Ritz vectors */
      /* Use a simple biorthogonalization that uses all vectors */   

      nAdded = nev;  //for the moment, we add all the vectors
      biortho_global_C(evecsl,ldvl,evecsr,ldvr,n,(*ncurEvals)+1,(*ncurEvals)+nev,3,params);
         
      //check the biorthogonality of the vectors
      /*
      for(i=0; i<(*ncurEvals)+nev; i++)
         for(k=0; k<(*ncurEvals)+nev; k++)
             {
                 tempc = wrap_zsum_cdot(&n,&evecsl[i*ldvl],&ONE,&evecsr[k*ldvr],&ONE,params);
                 fprintf(outputFile,"evecsl[%d]'*evecsr[%d]=%g %g\n",i,k,tempc.r,tempc.i);
             } 
      */

      /* Augument H */   
      /* (1:ncurEvals+nAdded,ncurEvals+1:ncurEvals+nAdded) block */
      for(k=(*ncurEvals); k<(*ncurEvals)+nAdded; k++)
          {
             matvec(&evecsr[k*ldvr],VR,params);
             for(i=0; i<(*ncurEvals)+nAdded; i++)
                {
                      tempc=wrap_zsum_cdot(&n,&evecsl[i*ldvl],&ONE,VR,&ONE,params);
                      H[i+k*ldh]=tempc;           
                }   
                   
          }


      /* (ncurEvals+1:ncurEvals+nAdded,1:ncurEvals) block */
      for(k=(*ncurEvals); k<(*ncurEvals)+nAdded; k++)
         {
             mathvec(&evecsl[k*ldvl],VL,params);
             for(i=0; i<(*ncurEvals); i++)
                {
                   tempc = wrap_zsum_cdot(&n,&evecsr[i*ldvr],&ONE,VL,&ONE,params);
                   H[k+i*ldh].r = tempc.r; 
                   H[k+i*ldh].i = -tempc.i;
                }
         } 
  
      (*ncurEvals) = (*ncurEvals) + nAdded;
   
        
      /* Reporting */
      wt2 = primme_get_wtime();
      if (plvl) {
            fprintf(outputFile, "Update\n");
            fprintf(outputFile, "Added %d vecs\n",nAdded);
            fprintf(outputFile, "U Wallclock : %-f\n", wt2-wt1);}
      
   } /* if phase==1 */

   /****************************************************/



   if(phase==2) //solve deflated bicgstab the correction equation
     {
        
        fprintf(outputFile,"\n\nDeflated bicgstab\n");;
        LRD_BICGSTAB_C(evecsl, ldvl,evecsr,ldvr,(*ncurEvals),&X[j*lde],b,     
                     lde,n,H,ldh,IPIV,work,matvec,params,(*AnormEst),maxit,(*restartTol),tol,ConvTestOpt,outputFile,&info);                
     } /*end of if(phase==2)*/
    
  } /*for(j=0; j<nrhs; j++) */





  // compute final evecs,etc.
  ComputeFinalEvecs_C
    ( (*ncurEvals),n,evecsl,ldvl,evecsr,ldvr,H,ldh,SRT_OPT,epsi, 
       evals, ernorms, xlnorms, xrnorms, angles, 
       matvec,params,work,6*lde);


  fprintf(outputFile,"\n\n Final Evals\n");
  for(i=0; i< (*ncurEvals); i++){
     fprintf(outputFile,"EVAL %-16.12E %-16.12E, ERNORM %-16.12E\n",
                         evals[i].r,evals[i].i,ernorms[i]);
     fprintf(outputFile,"XLNORM %-16.12E, XRNORM %-16.12E, ANGLE %-16.12E %-16.12E\n\n",
                         xlnorms[i],xrnorms[i],angles[i].r,angles[i].i);}
  fprintf(outputFile,"===================================================\n");

      
   return;
}
/******************************************************************************/
/* END OF IncrEigbicg_debug */
/******************************************************************************/

void init_BICG_C(Complex_C *vl, int ldvl, Complex_C *vr, int ldvr, int nvecs, 
               Complex_C *xinit, Complex_C *b, int lde, int n,
               Complex_C *H, int ldH, int *IPIV, Complex_C *work,
               void (*matvec) (void *, void *, void *), void *params, int *info)
{

    /* xinit = xinit + (vr inv(H) vl' res) */
    
    int i,ONE=1;
    Complex_C tempc,tpone={1.00000e+00,0.000000e+00};
    char cN='N';

    //compute res=b-Ax and store it in work[nvecs:nvecs+n-1]
    tempc = wrap_zsum_cdot(&n,xinit,&ONE,xinit,&ONE,params);
    if(tempc.r > 0.0) //if nonzero initial guess
      {
           matvec(xinit,&work[nvecs],params);

           for(i=0; i<n; i++)
              { work[nvecs+i].r = b[i].r - work[nvecs+i].r;
                work[nvecs+i].i = b[i].i - work[nvecs+i].i;}
      }
     else
       BLAS_CCOPY(&n,b,&ONE,&work[nvecs],&ONE);

    for(i=0; i<nvecs; i++){
      work[i] = wrap_zsum_cdot(&n, &vl[i*ldvl],&ONE,&work[nvecs],&ONE,params);}


    //solve using LU factorized H
    BLAS_CGETRS(&cN,&nvecs,&ONE,H,&ldH,IPIV,work,&nvecs,info);
    printf("inside init_bicg\n");

    if( (*info) != 0)
      { fprintf(stderr,"Error in BLAS_CGESV inside init-BICG. info %d\n",(*info));
        exit(2);
      }

    BLAS_CGEMV(&cN,&n,&nvecs,&tpone,vr,&ldvr,work,&ONE,&tpone,xinit,&ONE);

    return;
} 

/************************************************************************************/
/* END OF init_BICG */      
/************************************************************************************/

void LRD_BICGSTAB_C(Complex_C *vl, int ldvl, Complex_C *vr, int ldvr, int nvecs, Complex_C *x, Complex_C *b,     
               int lde, int n, Complex_C *H, int ldH, int *IPIV, Complex_C *work, void (*matvec) (void *, void *, void *), 
               void *params, float AnormEst, int maxiter, float DefTol, float tol, int ConvTestOpt, FILE *outputFile,  int *info)
{
     Complex_C *xincr, *resid, *tmpH;  //used to solve the correction equation
     int i,j,k,numIts,iters_used,cs,is,fs,allelems,ONE=1;
     Complex_C tzero={00.00e+00,00.00e+00};
     float resNorm,bnorm,curTol,leftoverTol;
     Complex_C tempc;
     float *reshist;   //residual norm history
     int flag1,flag2;
     float stoptol_used,stoptol_cur,xnorm,rhsnorm; 
     float MACHEPS=1.e-7;

     
     cs = sizeof(Complex_C);
     is = sizeof(int);
     fs = sizeof(float);

     xincr = (Complex_C *) calloc(lde,cs);
     if(xincr==NULL){
        printf("Error in allocating xincr in LRD_BICGSTAB\n");
        exit(1);
     }


     resid = (Complex_C *) calloc(lde,cs);
     if(resid==NULL){
        printf("Errro in allocating resid in LRD_BICGSTAB\n");
        exit(1);
     }




     tmpH = (Complex_C  *) calloc(ldH*nvecs,cs);  //temporary copy of H
     if(tmpH == NULL){
        printf("Error in allocating tmpH\n");
        exit(1);
     }

     reshist = (float *) calloc(maxiter,fs);
     if(reshist == NULL){
        printf("Error allocating reshist\n");
        exit(1);
     }

     
     //initial residual std::vector
     tempc = wrap_zsum_cdot(&n,x,&ONE,x,&ONE,params);
     xnorm = sqrt(tempc.r);
     if(tempc.r > 0 )
       {
          matvec(x,resid,params);
         for(i=0; i<n; i++){
             resid[i].r = b[i].r - resid[i].r;
             resid[i].i = b[i].i - resid[i].i;}
       }
     else
       BLAS_CCOPY(&n,b,&ONE,resid,&ONE);


     tempc = wrap_zsum_cdot(&n,resid,&ONE,resid,&ONE,params);
     resNorm = sqrt(tempc.r);

     tempc = wrap_zsum_cdot(&n,b,&ONE,b,&ONE,params);
     bnorm = sqrt(tempc.r);


     numIts=0; 
     //leftoverTol = DefTol;

     stoptol_used = tol*bnorm;
     if(ConvTestOpt==2)
        {
           stoptol_cur = MACHEPS*(AnormEst*xnorm + bnorm);
           if(stoptol_used < stoptol_cur)
               stoptol_used = stoptol_cur;
        }


     //Factorize the nvecs*nvecs part of H
     allelems = ldH*nvecs; 
     BLAS_CCOPY(&allelems,H,&ONE,tmpH,&ONE);
     BLAS_CGETRF(&nvecs,&nvecs,tmpH,&ldH,IPIV,info);
     
     if(info==0)
        {
           printf("ERROR: factorization of tmpH\n");
           exit(2);
        }

     if(DefTol < tol)
       curTol = tol;
     else
       curTol = DefTol;
     
     while( resNorm > stoptol_used){
     //int jrest=0;
     //while( jrest< 2){
     //    fprintf(outputFile,"restart\n");
         iters_used =0;         
         init_BICG_C(vl,ldvl,vr,ldvr,nvecs,xincr,resid,lde,n,
                   tmpH,ldH,IPIV,work,matvec,params,info);
         
         if( (*info) != 0){
            printf("Error in init_BICG\n");
            exit(2);
         } 


         
         fprintf(outputFile,"bnorm %g, curTol %g\n",bnorm,curTol);

         MyBICGSTAB_C(n,lde,xincr,resid,curTol,maxiter,&iters_used,reshist,
                     matvec,params,AnormEst, ConvTestOpt, work,&flag1);
         
         //if(flag1==0 || flag1=3)
         //   {
                 for(i=0; i<iters_used; i++){
                     fprintf(outputFile,"%-6d  %-22.16g \n",numIts+i,reshist[i]);} 

                 numIts = numIts + iters_used;

                 for(i=0; i<n; i++){
                     x[i].r = x[i].r + xincr[i].r;
                     x[i].i = x[i].i + xincr[i].i;
                     xincr[i]=tzero;}
                 

                 matvec(x,resid,params);
          
                 for(i=0; i<n; i++){
                     resid[i].r = b[i].r - resid[i].r;
                     resid[i].i = b[i].i - resid[i].i;}

                 tempc = wrap_zsum_cdot(&n, resid, &ONE, resid, &ONE, params);
                 resNorm = sqrt(tempc.r);
                 
                 //    }
         // else
             fprintf(outputFile,"BICGSTAB returns with flag %d\n",flag1); 
          
          if(ConvTestOpt ==2)
            {  
                 tempc = wrap_zsum_cdot(&n,x,&ONE,x,&ONE,params);
                 xnorm = sqrt(tempc.r);
                 stoptol_cur = MACHEPS*(AnormEst*xnorm+bnorm);
                 if(stoptol_used < stoptol_cur)
                    stoptol_used = stoptol_cur;
            }
         
          leftoverTol = stoptol_used/resNorm;
          if(leftoverTol < DefTol)
            curTol=DefTol;
          else
            curTol=leftoverTol;
          //fprintf(outputFile,"leftoverTol %g\n",leftoverTol);
          //jrest = jrest +1;
     }
     
     


     return;
}

/************************************************************************************************************/
                  
void ComputeFinalEvecs_C
    ( int nvecs, int n, Complex_C *evecsl, int ldvl, Complex_C *evecsr, int ldvr, Complex_C *H, int ldh, char SRT_OPT, 
      float epsi, Complex_C *evals, float *ernorms, float *xlnorms, float *xrnorms, Complex_C *angles, 
      void (*matvec)(void *, void *, void *), void *params, Complex_C *work, int worksize)
{

     Complex_C *tmpH,*COEFL,*COEFR,*Res;
     int i,allelems,infoCG,ONE=1;
     

     tmpH = (Complex_C *) calloc(ldh*nvecs,sizeof(Complex_C));
     if(tmpH == NULL){
        printf("ERROR: not able to allocate tmpH in ComputeFinalEvecs\n");
        exit(2);}


     COEFL = (Complex_C *) calloc(nvecs*nvecs,sizeof(Complex_C));
     if(COEFL == NULL){
        printf("ERROR: not able to allocate COEFL in ComputeFinalEvecs\n");
        exit(2);}

     COEFR = (Complex_C *) calloc(nvecs*nvecs,sizeof(Complex_C));
     if(COEFR == NULL){
        printf("ERROR: not able to allocate COEFR in ComputeFinalEvecs\n");
        exit(2);} 

     Res = (Complex_C *) calloc(n,sizeof(Complex_C));
     if(Res == NULL){
         printf("ERROR: not able to allocate Res in ComputeFinalEvecs\n");
         exit(2);}


     allelems=ldh*nvecs;
     BLAS_CCOPY(&allelems,H,&ONE,tmpH,&ONE);

     //void CG_eval(Complex_C *A, int N, int LDA, Complex_C *W, char SRT_OPT, 
     //        Complex_C *VL, int LDVL, Complex_C *VR, int LDVR, int *info);

     CG_eval(tmpH,nvecs,ldh,evals,SRT_OPT,epsi,COEFL,nvecs,COEFR,nvecs,&infoCG);

     if(infoCG != 0){
         printf("ERROR: CG_eval inside ComputeFinalEvecs returned with flag %d\n",infoCG);
         exit(3);}
        
      
     //void restart_X(Complex_C *X, int ldx, Complex_C *hVecs, int nLocal, 
     //          int basisSize, int restartSize, Complex_C *rwork, int rworkSize) 


     Crestart_X(evecsl, ldvl, COEFL, n, nvecs, nvecs, work, worksize);
     Crestart_X(evecsr, ldvr, COEFR, n, nvecs, nvecs, work, worksize); 


    //No need to biorthogonalize because evecsl and evecsr were orginally biorthogonal 
    //and COEFL and COEFR comes out of CG_eval as biorthogonal.

    //Compute Ritz values, residual norms, etc.
    for(i=0; i < nvecs; i++){
       //void computeResNorm( Complex_C *xr, Complex_C *xl, Complex_C *lambda, float *rnorm, int n, Complex_C *Res, 
       //float *xlnorm, float *xrnorm, Complex_C *cangle, void (*matvec)(void *, void *, void *), void *params)
       CcomputeResNorm(&evecsr[i*ldvr], &evecsl[i*ldvl],&evals[i],&ernorms[i],n,Res,&xlnorms[i],&xrnorms[i],
                      &angles[i],matvec,params);}

 

return;

}


/*********************DOUBLE PRECISION VERSION******************************************************/

void IncrEigbicg_Z(  int n, int lde,int nrhs, Complex_Z *X, Complex_Z *B, int *ncurEvals,         	
                     int ldh, Complex_Z *evecsl, Complex_Z *evecsr, Complex_Z *evals, 		
	             Complex_Z *H, void (*matvec) (void *, void *, void *),  
	             void (*mathvec)(void *, void *, void *), void *params, double *AnormEst, 
	             Complex_Z *work, Complex_Z *VL, int ldvl, Complex_Z *VR, int ldvr,        
	             Complex_Z *ework, int esize, double tol, double *restartTol, 	
	             int maxit, char SRT_OPT, double epsi, int ConvTestOpt, int plvl, int nev,
                     int v_max,FILE *outputFile)      
{

  /* Timing vars */
  double wt1,wt2,wE,wI;

  /* Pointers */
  Complex_Z  tempc, tempc1, tempc2, *tmpH, *x, *resid, *b;
  double     *rnorms, *reshist, normb, curTol,resNorm,leftTol;
  int        i,j,k, *IPIV, ONE = 1;
  int        zs, ds, tmpsize, is, phase, allelems ;
  int        numIts, flag, flag2,nAdded, nev_used, iters_used, info;
  Complex_Z  tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
  char       cR = 'R'; char cL = 'L'; char cN ='N'; 
  char       cV = 'V'; char cU = 'U'; char cC ='C';
  double     *xrnorms,*xlnorms,*ernorms;
  Complex_Z  *angles;
  

  /* ------------------------------------------------------------------- */
  /* Work allocations */
  /* ------------------------------------------------------------------- */ 
  zs = sizeof(Complex_Z); 
  ds = sizeof(double);
  is = sizeof(int);


  if( (IPIV = (int *) calloc(ldh,is)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate IPIV\n");
      exit(1);}


  if( (x = (Complex_Z *) calloc(lde,zs)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate x\n");
      exit(1);}


  if( (resid = (Complex_Z *) calloc(lde,zs)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate resid\n");
      exit(1);}

  
  if( (tmpH = (Complex_Z *) calloc(ldh*ldh,zs)) == NULL)
    { fprintf(stderr,"ERROR IncrEigbicg could not allocate tmpH\n");
      exit(1);}


  if ((work = (Complex_Z *) calloc(6*lde, zs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate work\n");
     exit(1);}
     
  if ((ework = (Complex_Z *) calloc(esize, zs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate ework\n");
     exit(1);}
  
  if ((VL = (Complex_Z *) calloc(v_max*ldvl, zs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate VL\n");
     exit(1);}

  if ((VR = (Complex_Z *) calloc(v_max*ldvr, zs)) == NULL) 
    {fprintf(stderr, "ERROR IncrEigbicg could not allocate VR\n");
     exit(1);}
     
  if ( (rnorms = (double *) calloc(ldh, ds))  == NULL )
     {fprintf(stderr, "ERROR IncrEigbicg could not allocate rnorms\n");
      exit(1);}

  if ( (reshist = (double *) calloc(maxit, ds)) == NULL)
     {fprintf(stderr, "ERROR IncrEigbicg could not allocate reshist\n");
      exit(1);}



  if( (xlnorms = (double *) calloc(ldh,ds)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate xlnorms\n");
     exit(1);}

  if( (xrnorms = (double *) calloc(ldh,ds)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate xrnorms\n");
     exit(1);}

  if( (ernorms = (double *) calloc(ldh,ds)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate ernorms\n");
     exit(1);}

  if( (angles = (Complex_Z *) calloc(ldh,zs)) == NULL){
     fprintf(stderr,"ERROR: IncrEigbicg couldn't allocate angles\n");
     exit(1);}

  /* ------------------------------------------------------------------- */
  /* end Work allocations */
  /* ------------------------------------------------------------------- */

  /* ---------------------------------------------------------------------------------  */
  /* Solving one by one the nrhs systems with incremental init-eigbicg or init-bicgstab */
  /* ---------------------------------------------------------------------------------- */


  for (j=0; j<nrhs; j++) {

      b = &B[j*lde];
      tempc = wrap_zdot(&n,b,&ONE,b,&ONE,params);
      printf("bnorm=%g\n",sqrt(tempc.r));
      normb = sqrt(tempc.r);

      numIts = 0;
      //choose eigbicg or bicgstab
      if(ldh-(*ncurEvals) >= nev )
        phase=1;
      else
        phase=2;

      if (plvl) fprintf(outputFile, "\n\nSystem %d\n", j);

      wE = 0.0; wI = 0.0;     /* Start accumulator timers */

      if ( (*ncurEvals > 0) && (phase==1)) {
         /* --------------------------------------------------------- */
         /* Perform init-BICG with evecsl and evecsr vectors          */
         /* xinit = xinit + evecsr*inv(H)*evecl'*(b-Ax0) 		      */
         /* --------------------------------------------------------- */
         wt1 = primme_get_wtime(); 
   
        /* copy H into tmpH otherwise it will be changed by CGEV */
        allelems = ldh*(*ncurEvals);
        BLAS_ZCOPY(&allelems,H,&ONE,tmpH,&ONE);
        // LU factorization of tmpH
        BLAS_ZGETRF(ncurEvals,ncurEvals,tmpH,&ldh,IPIV,&info);


        for(i=0; i<n; i++)
            x[i]=tzero;

        matvec(&X[j*lde],resid,params);

        for(i=0; i<n; i++){
         resid[i].r = b[i].r - resid[i].r;
         resid[i].i = b[i].i - resid[i].i;}


        init_BICG_Z(evecsl,ldvl,evecsr,ldvr,(*ncurEvals), 
                    x,resid,lde,n,tmpH,ldh,IPIV,work,matvec,params,&info);

        if(phase==1){
           for(i=0; i<n; i++){
              X[j*lde+i].r = X[j*lde+i].r + x[i].r;
              X[j*lde+i].i = X[j*lde+i].i + x[i].i;}}
        

        wt2 = primme_get_wtime();
	wI = wI + wt2-wt1;
       }
       /* end of init-BICG with evecsl and evecsr vectors              */
       /* ------------------------------------------------------------ */

       if(phase == 1){	
         /* ------------------------------------------------------------ */
         /* Solve Ax = b with x initial guess using eigbicg and compute
            new nev eigenvectors                                         */
         /* ------------------------------------------------------------ */
         wt1 = primme_get_wtime(); 
        
         Zeigbicg(n, lde, &X[j*lde], b, &normb, tol, maxit, SRT_OPT,epsi,ConvTestOpt, &numIts, reshist, 
                  &flag, plvl, work, matvec, mathvec, params, AnormEst, nev, 
                  &evals[(*ncurEvals)], &rnorms[(*ncurEvals)],
  	          v_max, VR,ldvr,VL,ldvl,esize,ework);

         wt2 = primme_get_wtime();
	 wE = wE + wt2-wt1;

         /* ---------- */
         /* Reporting  */
         /* ---------- */
         tempc1 = wrap_zdot(&n,&X[j*lde],&ONE,&X[j*lde],&ONE,params);
         if (plvl) {
            fprintf(outputFile, "For this rhs:\n");
            fprintf(outputFile, "Norm(solution) %-16.12E, AnormEst %-16.12E\n",sqrt(tempc1.r),(*AnormEst));
            fprintf(outputFile, "Total initBICG Wallclock : %-f\n", wI);
            fprintf(outputFile, "Total eigbicg Wallclock : %-f\n", wE);
            fprintf(outputFile, "Iterations: %-d\n", numIts); 
            fprintf(outputFile, "Actual Resid of LinSys  : %e\n", reshist[numIts-1]);
	    if (plvl > 1) 
               for (i=0; i < nev; i++) 
                   fprintf(outputFile, "Eval[%d]: %-22.15E    %-22.15E           rnorm: %-22.15E\n", 
	                       i+1, evals[*ncurEvals+i].r,  evals[*ncurEvals+i].i, rnorms[*ncurEvals+i]); 
            if (plvl >1)
              {
                 fprintf(outputFile,"Residual norm\n");
               for( i=0; i < numIts; i++)
                   fprintf(outputFile,"%-d  %-22.15E\n",i,reshist[i]);
              }  

            if (flag != 0) {
               fprintf(outputFile, "Error: eigbicg returned with nonzero exit status\n");
            return;}
         }

      /* ------------------------------------------------------------------- */
      /* Update the evecsl, evecsr,  and evecsl'*A*evecsr                               */
      /* ------------------------------------------------------------------- */
      wt1 = primme_get_wtime(); 


      /* Append new Ritz pairs to evecs */
      for (i=0; i<nev; i++){
	 BLAS_ZCOPY(&n, &VL[i*ldvl], &ONE, &evecsl[((*ncurEvals)+i)*ldvl], &ONE);
	 BLAS_ZCOPY(&n, &VR[i*ldvr], &ONE, &evecsr[((*ncurEvals)+i)*ldvr], &ONE);}
  
      /* Bi-Orthogonalize the new Ritz vectors */
      /* Use a simple biorthogonalization that uses all vectors */   

      nAdded = nev;  //for the moment, we add all the vectors
      biortho_global_Z(evecsl,ldvl,evecsr,ldvr,n,(*ncurEvals)+1,(*ncurEvals)+nev,3,params);
         
      //check the biorthogonality of the vectors
      /*
      for(i=0; i<(*ncurEvals)+nev; i++)
         for(k=0; k<(*ncurEvals)+nev; k++)
             {
                 tempc = wrap_zdot(&n,&evecsl[i*ldvl],&ONE,&evecsr[k*ldvr],&ONE,params);
                 fprintf(outputFile,"evecsl[%d]'*evecsr[%d]=%g %g\n",i,k,tempc.r,tempc.i);
             } 
      */

      /* Augument H */   
      /* (1:ncurEvals+nAdded,ncurEvals+1:ncurEvals+nAdded) block */
      for(k=(*ncurEvals); k<(*ncurEvals)+nAdded; k++)
          {
             matvec(&evecsr[k*ldvr],VR,params);
             for(i=0; i<(*ncurEvals)+nAdded; i++)
                {
                      tempc=wrap_zdot(&n,&evecsl[i*ldvl],&ONE,VR,&ONE,params);
                      H[i+k*ldh]=tempc;           
                }   
                   
          }


      /* (ncurEvals+1:ncurEvals+nAdded,1:ncurEvals) block */
      for(k=(*ncurEvals); k<(*ncurEvals)+nAdded; k++)
         {
             mathvec(&evecsl[k*ldvl],VL,params);
             for(i=0; i<(*ncurEvals); i++)
                {
                   tempc = wrap_zdot(&n,&evecsr[i*ldvr],&ONE,VL,&ONE,params);
                   H[k+i*ldh].r = tempc.r; 
                   H[k+i*ldh].i = -tempc.i;
                }
         } 
  
      (*ncurEvals) = (*ncurEvals) + nAdded;
   
        
      /* Reporting */
      wt2 = primme_get_wtime();
      if (plvl) {
            fprintf(outputFile, "Update\n");
            fprintf(outputFile, "Added %d vecs\n",nAdded);
            fprintf(outputFile, "U Wallclock : %-f\n", wt2-wt1);}
      
   } /* if phase==1 */

   /****************************************************/



   if(phase==2) //solve deflated bicgstab the correction equation
     {
        
        fprintf(outputFile,"\n\nDeflated bicgstab\n");;
        LRD_BICGSTAB_Z(evecsl, ldvl,evecsr,ldvr,(*ncurEvals),&X[j*lde],b,     
                     lde,n,H,ldh,IPIV,work,matvec,params,(*AnormEst),maxit,(*restartTol),tol,ConvTestOpt,outputFile,&info);                
     } /*end of if(phase==2)*/
    
  } /*for(j=0; j<nrhs; j++) */





  // compute final evecs,etc.
  ComputeFinalEvecs_Z
    ( (*ncurEvals),n,evecsl,ldvl,evecsr,ldvr,H,ldh,SRT_OPT,epsi, 
       evals, ernorms, xlnorms, xrnorms, angles, 
       matvec,params,work,6*lde);


  fprintf(outputFile,"\n\n Final Evals\n");
  for(i=0; i< (*ncurEvals); i++){
     fprintf(outputFile,"EVAL %-16.12E %-16.12E, ERNORM %-16.12E\n",
                         evals[i].r,evals[i].i,ernorms[i]);
     fprintf(outputFile,"XLNORM %-16.12E, XRNORM %-16.12E, ANGLE %-16.12E %-16.12E\n\n",
                         xlnorms[i],xrnorms[i],angles[i].r,angles[i].i);}
  fprintf(outputFile,"===================================================\n");

      
   return;
}
/******************************************************************************/
/* END OF IncrEigbicg_Z */
/******************************************************************************/

void init_BICG_Z(Complex_Z *vl, int ldvl, Complex_Z *vr, int ldvr, int nvecs, 
               Complex_Z *xinit, Complex_Z *b, int lde, int n,
               Complex_Z *H, int ldH, int *IPIV, Complex_Z *work,
               void (*matvec) (void *, void *, void *), void *params, int *info)
{

    /* xinit = xinit + (vr inv(H) vl' res) */
    
    int i,ONE=1;
    Complex_Z tempc,tpone={1.00000e+00,0.000000e+00};
    char cN='N';

    //compute res=b-Ax and store it in work[nvecs:nvecs+n-1]
    tempc = wrap_zdot(&n,xinit,&ONE,xinit,&ONE,params);
    if(tempc.r > 0.0) //if nonzero initial guess
      {
           matvec(xinit,&work[nvecs],params);

           for(i=0; i<n; i++)
              { work[nvecs+i].r = b[i].r - work[nvecs+i].r;
                work[nvecs+i].i = b[i].i - work[nvecs+i].i;}
      }
     else
       BLAS_ZCOPY(&n,b,&ONE,&work[nvecs],&ONE);

    for(i=0; i<nvecs; i++){
      work[i] = wrap_zdot(&n, &vl[i*ldvl],&ONE,&work[nvecs],&ONE,params);}


    //solve using LU factorized H
    BLAS_ZGETRS(&cN,&nvecs,&ONE,H,&ldH,IPIV,work,&nvecs,info);
    printf("inside init_bicg\n");

    if( (*info) != 0)
      { fprintf(stderr,"Error in BLAS_ZGESV inside init-BICG_Z. info %d\n",(*info));
        exit(2);
      }

    BLAS_ZGEMV(&cN,&n,&nvecs,&tpone,vr,&ldvr,work,&ONE,&tpone,xinit,&ONE);

    return;
} 

/************************************************************************************/
/* END OF init_BICG      
************************************************************************************/

void LRD_BICGSTAB_Z(Complex_Z *vl, int ldvl, Complex_Z *vr, int ldvr, int nvecs, Complex_Z *x, Complex_Z *b,     
                    int lde, int n, Complex_Z *H, int ldH, int *IPIV, Complex_Z *work, void (*matvec) (void *, void *, void *), 
                    void *params, double AnormEst, int maxiter, double DefTol, double tol, int ConvTestOpt, FILE *outputFile,  int *info)
{
     Complex_Z *xincr, *resid, *tmpH;  //used to solve the correction equation
     int i,j,k,numIts,iters_used,zs,is,ds,allelems,ONE=1;
     Complex_Z tzero={00.00e+00,00.00e+00};
     double resNorm,bnorm,curTol,leftoverTol;
     Complex_Z tempc;
     double *reshist;   //residual norm history
     int flag1,flag2;
     double stoptol_used,stoptol_cur,xnorm,rhsnorm; 
     double MACHEPS=1.e-16;

     
     zs = sizeof(Complex_Z);
     is = sizeof(int);
     ds = sizeof(double);

     xincr = (Complex_Z *) calloc(lde,zs);
     if(xincr==NULL){
        printf("Error in allocating xincr in LRD_BICGSTAB\n");
        exit(1);
     }


     resid = (Complex_Z *) calloc(lde,zs);
     if(resid==NULL){
        printf("Errro in allocating resid in LRD_BICGSTAB\n");
        exit(1);
     }




     tmpH = (Complex_Z  *) calloc(ldH*nvecs,zs);  //temporary copy of H
     if(tmpH == NULL){
        printf("Error in allocating tmpH\n");
        exit(1);
     }

     reshist = (double *) calloc(maxiter,ds);
     if(reshist == NULL){
        printf("Error allocating reshist\n");
        exit(1);
     }

     
     //initial residual std::vector
     tempc = wrap_zdot(&n,x,&ONE,x,&ONE,params);
     xnorm = sqrt(tempc.r);
     if(tempc.r > 0 )
       {
          matvec(x,resid,params);
         for(i=0; i<n; i++){
             resid[i].r = b[i].r - resid[i].r;
             resid[i].i = b[i].i - resid[i].i;}
       }
     else
       BLAS_ZCOPY(&n,b,&ONE,resid,&ONE);


     tempc = wrap_zdot(&n,resid,&ONE,resid,&ONE,params);
     resNorm = sqrt(tempc.r);

     tempc = wrap_zdot(&n,b,&ONE,b,&ONE,params);
     bnorm = sqrt(tempc.r);


     numIts=0; 
     //leftoverTol = DefTol;

     stoptol_used = tol*bnorm;
     if(ConvTestOpt==2)
       {
           stoptol_cur = MACHEPS*(AnormEst*xnorm + bnorm);
           if(stoptol_used < stoptol_cur)
               stoptol_used = stoptol_cur;
       }


     //Factorize the nvecs*nvecs part of H
     allelems = ldH*nvecs; 
     BLAS_ZCOPY(&allelems,H,&ONE,tmpH,&ONE);
     BLAS_ZGETRF(&nvecs,&nvecs,tmpH,&ldH,IPIV,info);
     
     if(info==0)
        {
           printf("ERROR: factorization of tmpH\n");
           exit(2);
        }

     if(DefTol < tol)
       curTol = tol;
     else
       curTol = DefTol;
     
     while( resNorm > stoptol_used){
     //int jrest=0;
     //while( jrest< 2){
     //    fprintf(outputFile,"restart\n");
         iters_used =0;         
         init_BICG_Z(vl,ldvl,vr,ldvr,nvecs,xincr,resid,lde,n,
                     tmpH,ldH,IPIV,work,matvec,params,info);
         
         if( (*info) != 0){
            printf("Error in init_BICG\n");
            exit(2);
         } 


         
         fprintf(outputFile,"bnorm %g, curTol %g\n",bnorm,curTol);

         MyBICGSTAB_Z(n,lde,xincr,resid,curTol,maxiter,&iters_used,reshist,
                     matvec,params,AnormEst, ConvTestOpt, work,&flag1);
         
         //if(flag1==0 || flag1=3)
         //   {
                 for(i=0; i<iters_used; i++){
                     fprintf(outputFile,"%-6d  %-22.16g \n",numIts+i,reshist[i]);} 

                 numIts = numIts + iters_used;

                 for(i=0; i<n; i++){
                     x[i].r = x[i].r + xincr[i].r;
                     x[i].i = x[i].i + xincr[i].i;
                     xincr[i]=tzero;}
                 

                 matvec(x,resid,params);
          
                 for(i=0; i<n; i++){
                     resid[i].r = b[i].r - resid[i].r;
                     resid[i].i = b[i].i - resid[i].i;}

                 tempc = wrap_zdot(&n, resid, &ONE, resid, &ONE, params);
                 resNorm = sqrt(tempc.r);
                 
                     //    }
         // else
             fprintf(outputFile,"BICGSTAB returns with flag %d\n",flag1); 

          if(ConvTestOpt==2)
            { 
                tempc = wrap_zdot(&n,x,&ONE,x,&ONE,params);
                xnorm = sqrt(tempc.r);
                stoptol_cur = MACHEPS*(AnormEst*xnorm+bnorm);
                if(stoptol_used < stoptol_cur)
                    stoptol_used = stoptol_cur;
             }
          leftoverTol = stoptol_used/resNorm;
          if(leftoverTol < DefTol)
            curTol=DefTol;
          else
            curTol=leftoverTol;
          fprintf(outputFile,"leftoverTol %g\n",leftoverTol);
          //jrest = jrest +1;
     }
     
     


     return;
}

/************************************************************************************************************/
                  
void ComputeFinalEvecs_Z
    ( int nvecs, int n, Complex_Z *evecsl, int ldvl, Complex_Z *evecsr, int ldvr, Complex_Z *H, int ldh, char SRT_OPT, 
      double epsi, Complex_Z *evals, double *ernorms, double *xlnorms, double *xrnorms, Complex_Z *angles, 
      void (*matvec)(void *, void *, void *), void *params, Complex_Z *work, int worksize)
{

     Complex_Z *tmpH,*COEFL,*COEFR,*Res;
     int i,allelems,infoCG,ONE=1;
     

     tmpH = (Complex_Z *) calloc(ldh*nvecs,sizeof(Complex_Z));
     if(tmpH == NULL){
        printf("ERROR: not able to allocate tmpH in ComputeFinalEvecs\n");
        exit(2);}


     COEFL = (Complex_Z *) calloc(nvecs*nvecs,sizeof(Complex_Z));
     if(COEFL == NULL){
        printf("ERROR: not able to allocate COEFL in ComputeFinalEvecs\n");
        exit(2);}

     COEFR = (Complex_Z *) calloc(nvecs*nvecs,sizeof(Complex_Z));
     if(COEFR == NULL){
        printf("ERROR: not able to allocate COEFR in ComputeFinalEvecs\n");
        exit(2);} 

     Res = (Complex_Z *) calloc(n,sizeof(Complex_Z));
     if(Res == NULL){
         printf("ERROR: not able to allocate Res in ComputeFinalEvecs\n");
         exit(2);}


     allelems=ldh*nvecs;
     BLAS_ZCOPY(&allelems,H,&ONE,tmpH,&ONE);

     //void CG_eval(Complex_C *A, int N, int LDA, Complex_C *W, char SRT_OPT, 
     //        Complex_C *VL, int LDVL, Complex_C *VR, int LDVR, int *info);

     ZG_eval(tmpH,nvecs,ldh,evals,SRT_OPT,epsi,COEFL,nvecs,COEFR,nvecs,&infoCG);

     if(infoCG != 0){
         printf("ERROR: CG_eval inside ComputeFinalEvecs returned with flag %d\n",infoCG);
         exit(3);}
        
      
     //void restart_X(Complex_C *X, int ldx, Complex_C *hVecs, int nLocal, 
     //          int basisSize, int restartSize, Complex_C *rwork, int rworkSize) 


     Zrestart_X(evecsl, ldvl, COEFL, n, nvecs, nvecs, work, worksize);
     Zrestart_X(evecsr, ldvr, COEFR, n, nvecs, nvecs, work, worksize); 


    //No need to biorthogonalize because evecsl and evecsr were orginally biorthogonal 
    //and COEFL and COEFR comes out of CG_eval as biorthogonal.

    //Compute Ritz values, residual norms, etc.
    for(i=0; i < nvecs; i++){
       //void computeResNorm( Complex_C *xr, Complex_C *xl, Complex_C *lambda, float *rnorm, int n, Complex_C *Res, 
       //float *xlnorm, float *xrnorm, Complex_C *cangle, void (*matvec)(void *, void *, void *), void *params)
       ZcomputeResNorm(&evecsr[i*ldvr], &evecsl[i*ldvl],&evals[i],&ernorms[i],n,Res,&xlnorms[i],&xrnorms[i],
                      &angles[i],matvec,params);}

 

return;

}

    
    


