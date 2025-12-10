/*******************************************************************************
    Eig-BiCG algorithm for solving ono-symmetric linear system Ax=b and 
    approximating few small eignvalues and eignvectors.
    
    Authors: Abdou M. Abdel-Rehim, Andreas Stathopolous , and Kostas Originos
             andreas@cs.wm.edu

    Last updated: 10/14/2009
*******************************************************************************/ 

/*-------------------------------------------------------------------------

   BLAS_<name> are macros that in numerical.h are defined to the
   appropriate way to call the Fortran BLAS <name> function. Typically 
   this involves an underscore at the end, eg. BLAS_saxpy -> saxpy_

   wrap_<name> is used for two functions. dot products and gemv for large
   parallel distributed vectors. The wrappers include a globalSum of the result

--------------------------------------------------------------------------*/


#include "qdp-lapack_eigbicg.h"


/* EIGBICG *************************************************************************************/
void Ceigbicg(int n, int lde, Complex_C *x, Complex_C *b, float *normb, float tol, int maxit, 
              char SRT_OPT, float epsi, int ConvTestOpt, int *iter, float *reshist, int *flag, 
              int plvl, Complex_C *work, void (*matvec) (void *, void *, void *), 
              void (*mathvec)(void *, void *, void *), void *params, float *AnormEst, int nev, 
              Complex_C *evals, float *rnorms, int v_max, Complex_C *VR, int LDVR, Complex_C *VL, 
              int LDVL, int esize, Complex_C *ework)
{

  
  Complex_C *r, *rs, *p, *ps, *Ap, *Atps, *resid;                     
  Complex_Z alpha, alphaprev, beta, betaprev, rho, rhoprev;           
  double    rnorm, rnormprev;                                 
  int       i, j,k,it, v_size,tmpi,allelems,info,ii,jj;
  int       zs, cs, ds, tmpsize;                              
  Complex_C tmpc,tmpc1,tmpc2,tmpc3,cof1,cof2;
  Complex_Z tmpz,tmpz1;  

  Complex_Z *H,*tmpH,*Hevals,*Hevecsl, *Hevecsr;            
  Complex_Z *Hevalsold,*Hevecslold, *Hevecsrold;     
  Complex_Z *ZCoefl, *ZCoefr;              
  Complex_C *Coefl , *Coefr;               

  Complex_C  *Ap_prev, *Atps_prev, *tmp_work;                   /* pointers in ework */
  float      tmpd,tmpd1,tmpd2,ernorm,xlnorm,xrnorm, xnorm, stoptol_used, stoptol_cur;
  Complex_C  lambda,cangle;


  float MACHEPS=1e-7;

  /* constants */
  int ONE = 1;      /* for passing 1 to BLAS routines */
  Complex_C tpone = {+1.00000000e+00,+0.0000000e00}, 
            tzero = {+0.00000000e+00,+0.0000000e00}; /*single precision Complex numbers 1 and zero */
  Complex_Z z_one = {+1.0000000000000000e+00,+0.0000000000000000e00}, 
            z_zero= {+0.0000000000000000e+00,+0.0000000000000000e00}; /*double precision Complex numbers 1 and zero */
  char cR = 'R'; char cL = 'L'; char cN ='N';                       /* Character constants */
  char cV = 'V'; char cU = 'U'; char cC ='C';

  //FILE *OFile=fopen("Ceigbicg_debug.dat","w");
  Complex_C tmpcr,tmpcl;
  float tmpfr,tmpfl;

  /* allocations */
  zs = sizeof(Complex_Z); 
  cs = sizeof(Complex_C); 
  ds = sizeof(double);


  if ((H = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
       fprintf(stderr, "ERROR Could not allocate H\n");
  if ((tmpH = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
       fprintf(stderr, "ERROR Could not allocate tmpH\n");
  if ((Hevecsl = (Complex_Z *)calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsl\n");
  if ((Hevecsr = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsr\n");
  if ((Hevecslold = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecslold\n");
  if ((Hevecsrold = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsrold\n");
  if ((Hevals = (Complex_Z *) calloc(v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevals\n");
  if ((Hevalsold = (Complex_Z *) calloc(v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevalsold\n");
  if ((Coefl = (Complex_C *) calloc(2*nev*v_max, cs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Coefl\n");
  if ((Coefr = (Complex_C *) calloc(2*nev*v_max, cs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Coefr\n");
  if ((ZCoefl = (Complex_Z *) calloc(2*nev*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate ZCoefl\n");
  if ((ZCoefr = (Complex_Z *) calloc(2*nev*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate ZCoefr\n");


  

  // Initialize H and tmpH to zeros
  for(i=0; i<v_max*v_max; i++)
     { H[i]=z_zero; tmpH[i]=z_zero;} 
 

  /* setup pointers into work */
  r    = work;
  rs   = work + lde;
  p    = work + 2*lde;
  ps   = work + 3*lde;
  Ap   = work + 4*lde;
  Atps = work + 5*lde;


  /* setup pointers into ework */
  Ap_prev  = ework;
  Atps_prev= ework+lde;
  tmp_work = ework+2*lde;
  tmpsize  = esize-2*lde;  /*leftover in tmp_work*/
  

  /*Initialization*/
  tmpc = wrap_zsum_cdot(&n, b, &ONE, b, &ONE, params); /* Norm of rhs, b */
  (*normb) = sqrt(tmpc.r);


  /* If right hand side is zero return zero solution. ITER=0 */
  if ((*normb) == 0.0) {		
     for (i=0; i<n; i++) 
             x[i]=tzero; 
     (*flag) = 0;
     (*iter) = 0; 		
     reshist[0] = 0.0;
     if (plvl) CdisplayInfo(tol,maxit,*flag,*iter,reshist[0]);
     return;
  }


  /* Zero-th residual: r = b - A*x  */
  matvec(x, r, params);			
  for (i = 0; i < n; i ++) {
      r[i].r = b[i].r - r[i].r;
      r[i].i = b[i].i - r[i].i;
   }

  /* Zero-th rs, p, and ps = Zero-th r */
  BLAS_CCOPY(&n, r, &ONE, rs, &ONE);
  BLAS_CCOPY(&n, r, &ONE, p, &ONE);
  BLAS_CCOPY(&n, rs, &ONE, ps, &ONE);

  tmpc = wrap_zsum_cdot(&n, r, &ONE, r, &ONE, params); /* Sum in double */
  reshist[0] = sqrt(tmpc.r); rnormprev=reshist[0];
  if(plvl) printf("Initial residual %g\n",reshist[0]); 

  tmpc = wrap_zsum_cdot(&n, rs, &ONE, r, &ONE, params); /* Sum in double */
  rhoprev.r =tmpc.r;  rhoprev.i=tmpc.i;
  rho=z_zero; alpha=z_zero; alphaprev=z_zero; beta=z_zero; betaprev=z_zero;



   
  v_size =0;   /* counts vectors in the eignstd::vector search spaces */
  it=0;        /* iteration counter*/

  while( it < maxit-1 ){      /*main bicg loop*/


      it = it + 1;
      matvec(p,Ap,params);
      mathvec(ps,Atps,params);

      if(v_size == v_max-1){
        BLAS_CCOPY(&n, Ap, &ONE, Ap_prev, &ONE);
        BLAS_CCOPY(&n, Atps, &ONE, Atps_prev, &ONE);}
      
      if( (v_size == v_max) && nev >0 ){  /*eigenvalue part*/
         
         allelems=v_max*v_max;
         BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE); //copy H into a temporary storage to be used in ZG_eval  

         ZG_eval(tmpH, v_max, v_max, Hevals, SRT_OPT, (double) epsi, 
                 Hevecsl, v_max, Hevecsr, v_max, &info);


         if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of H(v_max,v_max)\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            (*flag) = 3;
            return; }


        //Estimated norm(A)
        tmpd = z_abs_primme(Hevals[v_max-1]);
        if(tmpd > (*AnormEst) )
          (*AnormEst) = tmpd;
        

         BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE);
         tmpi=v_max-1;  

         ZG_eval(tmpH,tmpi,v_max,Hevalsold,SRT_OPT, (double) epsi, Hevecslold,v_max,
                 Hevecsrold,v_max,&info);



         if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of H(v_max-1,v_max-1)\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            *flag = 3;
            return; }
          

         /* set the v_max elements of Hevecslold and Hevecsrold to zero */
         for(i=1; i <= tmpi ; i++){
           Hevecslold[i*v_max-1].r = 0.0 ; Hevecslold[i*v_max-1].i = 0.0;
           Hevecsrold[i*v_max-1].r = 0.0 ; Hevecsrold[i*v_max-1].i = 0.0;}

        /* append the nev eigenvectors computed at step (v_max-1) after the nev
           eigenvectors computed at step v_max */


        tmpi=nev*v_max;
        BLAS_ZCOPY(&tmpi, Hevecsrold, &ONE, &Hevecsr[tmpi], &ONE);
        BLAS_ZCOPY(&tmpi, Hevecslold, &ONE, &Hevecsl[tmpi], &ONE);



        /* bi-orthogonalize the nev old vectors againist the new nev vectors */
        biortho_local(Hevecsl, v_max, Hevecsr, v_max, v_max, nev+1, 2*nev, 3);
   


        /*Compute new matrix H(1:2nev,1:2nev) = Hevecl(1:v_max,1:2nev)'*H*Hevecr(1:v_max,1:2nev) */
        /* first compute H*Hevecsr and store it it Hevecsrold */

        tmpi=2*nev;
        BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &v_max, &z_one, H, &v_max,
                    Hevecsr, &v_max, &z_zero, Hevecsrold, &v_max);   
        /* next compute Hevecsl(1:v_max,1:2nev)' * Hevecsrold(1:v_max, 1:2nev) 
           and store the result in tmpH such that the leading dimension
           of tmpH will be 2*nev and the actual result will be in the 
           2nev*2nev block of tmpH */
        BLAS_ZGEMM( &cC, &cN, &tmpi, &tmpi, &v_max, &z_one, Hevecsl,
                    &v_max, Hevecsrold, &v_max, &z_zero, tmpH, &tmpi);

       /* solve the 2nev*2nev eigenvalue problem */
       /* since tmpH won't be needed afterwards, we can pass it to GEEV 
          store eigenvalues and eigenvectors in the old arrays Hevalsold,.etc*/

       ZG_eval(tmpH,tmpi,tmpi,Hevalsold,SRT_OPT,(double) epsi, 
               Hevecslold,tmpi,Hevecsrold,tmpi,&info);


       /* check that the eigenvalue problem didn't fail */
       if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of the 2nevx2nev problem\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            *flag=3;
            return; }


       /* Compute the coefficient matrices */
       BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &tmpi, &z_one, Hevecsr,
                   &v_max, Hevecsrold, &tmpi, &z_zero, ZCoefr, &v_max);

       
       BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &tmpi, &z_one, Hevecsl,
                   &v_max, Hevecslold, &tmpi, &z_zero, ZCoefl, &v_max);





       /* Copy double precision Coef matrices into single precision ones */
       allelems=2*nev*v_max;
       for (i=0; i< allelems; i++){
          Coefr[i].r=ZCoefr[i].r;  Coefr[i].i=ZCoefr[i].i;
          Coefl[i].r=ZCoefl[i].r;  Coefl[i].i=ZCoefl[i].i;}



       /* Restarting */
       v_size=2*nev;

       Crestart_X(VL, LDVL, Coefl, n, v_max, v_size, tmp_work, tmpsize);
       Crestart_X(VR, LDVR, Coefr, n, v_max, v_size, tmp_work, tmpsize);

       //-----------DEBUG---------------------------
       //Nromalize the new 2*nev vectors VL and VR:
       /*
       for(i=0; i<v_size; i++)
         {
             tmpc = wrap_zsum_cdot(&n,&VR[i*LDVR],&ONE,&VR[i*LDVR],&ONE,params);
             for(j=0; j<n; j++)
                {
                    VR[j+i*LDVR].r = VR[j+i*LDVR].r / sqrt(tmpc.r);
                    VR[j+i*LDVR].i = VR[j+i*LDVR].i / sqrt(tmpc.r);
                }
             tmpc = wrap_zsum_cdot(&n,&VR[i*LDVR],&ONE,&VL[i*LDVL],&ONE,params);
             c_div_primme(&tmpc1,&tpone,&tmpc);             
             BLAS_CSCAL(&n,&tmpc1,&VL[i*LDVL],&ONE);
         }
      */
      //-------------END DEBUG-----------------------
     



       /* Restart H = diag(Hevals) plus a column and a row */
       allelems = v_max*v_max;
       for (i = 0; i < allelems; i++ )  
             H[i]= z_zero; 

       for (i = 0; i < v_size; i++) {
           H[i*(v_max+1)].r = Hevalsold[i].r;
           H[i*(v_max+1)].i = Hevalsold[i].i;}
           

       /* Compute the elemnts 2nev+1 row and column of H */
       tmpc.r = beta.r; tmpc.i=-beta.i ;  /* beta' to a single precision: WARNING */
       BLAS_CSCAL( &n, &tmpc, Atps_prev, &ONE);  /* beta'*Atps_prev */
       
       for(i=0; i<n; i++){
          Atps_prev[i].r = Atps[i].r - Atps_prev[i].r;
          Atps_prev[i].i = Atps[i].i - Atps_prev[i].i;}

       tmpc.r = beta.r; tmpc.i= beta.i;
       BLAS_CSCAL( &n, &tmpc, Ap_prev, &ONE);

       for(i=0; i<n; i++){
          Ap_prev[i].r = Ap[i].r - Ap_prev[i].r;
          Ap_prev[i].i = Ap[i].i - Ap_prev[i].i;}

       z_div_primme(&tmpz1,&z_one,&rho);   // 1/rho

  
       for(i=0; i<v_size; i++){
          tmpc= wrap_zsum_cdot(&n, Atps_prev, &ONE, &VR[i*LDVR], &ONE, params); 
          H[i*v_max+v_size].r= rnorm*(tmpc.r*tmpz1.r-tmpc.i*tmpz1.i); 
          H[i*v_max+v_size].i= rnorm*(tmpc.r*tmpz1.i+tmpc.i*tmpz1.r);
          tmpc = wrap_zsum_cdot(&n, &VL[i*LDVL], &ONE, Ap_prev, &ONE, params);
          H[i+v_size*v_max].r = tmpc.r/rnorm; H[i+v_size*v_max].i= tmpc.i/rnorm;}

       if(plvl>=4)
         {

             for(i=0; i<2*nev; i++)
               {
                  
                  //temprarily use Ap_prev to store the eigenvalue residual std::vector.
                  //Ap_prev is not needed at this point till next restart when it will be
                  //recomputed.
                  CcomputeResNorm(&VR[i*LDVR],&VL[i*LDVL],&lambda,&ernorm,n,Ap_prev, 
                                 &xlnorm,&xrnorm,&cangle,matvec,params);
                  printf("%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                          i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);

                  //fprintf(OFile,"%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                  //        i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);




               }
         } 


      } /*end of the eigenvalue part if(v_size==v_max & nev >0) */

      /* continue bicg from that point on */
      
      // add a new std::vector to the basis

      if(nev > 0)
        {     
           for(i=0; i<n; i++) {
              VL[LDVL*v_size + i].r = rs[i].r;  VL[LDVL*v_size+i].i = rs[i].i;
              VR[LDVR*v_size+i].r   = r[i].r ;  VR[LDVR*v_size+i].i = r[i].i ; } 
      
           tmpz.r=rnormprev; tmpz.i=0.0e+00;
           z_div_primme( &tmpz1, &z_one, &tmpz);
           tmpc.r = tmpz1.r;   tmpc.i = tmpz1.i;
           BLAS_CSCAL( &n, &tmpc, &VR[LDVR*v_size], &ONE);

           z_div_primme( &tmpz1, &tmpz, &rhoprev);
           tmpc.r= tmpz1.r;  tmpc.i = - tmpz1.i;
           BLAS_CSCAL( &n, &tmpc, &VL[LDVL*v_size], &ONE);
           v_size++;
        }

      tmpc = wrap_zsum_cdot(&n, ps, &ONE, Ap, &ONE, params);

      tmpz.r=tmpc.r; tmpz.i=tmpc.i;

      z_div_primme( &alpha, &rhoprev, &tmpz);  
      
      if(z_abs_primme(alpha) < 1.e-14 ) { *flag=2; break;}
      
      tmpc.r = -alpha.r; tmpc.i = -alpha.i;      //update r
      BLAS_CAXPY( &n, &tmpc, Ap, &ONE, r, &ONE);

      tmpc.r= -alpha.r; tmpc.i = alpha.i;
      BLAS_CAXPY( &n, &tmpc, Atps, &ONE, rs, &ONE);  //update rs

      tmpc = wrap_zsum_cdot(&n, r, &ONE, r, &ONE, params); /* Sum in double */
      rnorm = sqrt(tmpc.r);
      reshist[it] = rnorm;
      if(plvl>=2) printf("iter# %d, residual norm = %0.12g\n",it,rnorm);

      tmpc.r = alpha.r; tmpc.i = alpha.i;          //update the solution
      BLAS_CAXPY( &n, &tmpc, p, &ONE, x, &ONE);




      tmpc= wrap_zsum_cdot(&n, rs, &ONE, r, &ONE, params);
      rho.r = tmpc.r;  rho.i = tmpc.i;
      if(z_abs_primme(rho)< 1.e-40){ *flag=2; break;}
 
      z_div_primme( &beta, &rho, &rhoprev);
      if(z_abs_primme(beta)< 1.e-14){ *flag=2; break;}
      


      tmpc.r= beta.r; tmpc.i = beta.i;           //update p
      BLAS_CSCAL( &n, &tmpc, p, &ONE);
      BLAS_CAXPY( &n, &tpone, r, &ONE, p, &ONE);

      tmpc.r = beta.r; tmpc.i = -beta.i;         //update ps
      BLAS_CSCAL( &n, &tmpc, ps, &ONE);
      BLAS_CAXPY( &n, &tpone, rs, &ONE, ps, &ONE);


      if(nev > 0)
        {
            /* elements of the projection matrix */
            /* diagonal elemnts */
      
            z_div_primme( &tmpz, &z_one, &alpha);
            if((v_size-1) == 0){ 
                H[(v_max+1)*(v_size-1)].r = tmpz.r;
                H[(v_max+1)*(v_size-1)].i = tmpz.i;}
            else{
                z_div_primme( &tmpz1, &betaprev, &alphaprev);
                H[(v_max+1)*(v_size-1)].r =  tmpz.r + tmpz1.r;
                H[(v_max+1)*(v_size-1)].i =  tmpz.i + tmpz1.i;}

            /* off-diagonal */
            if( v_size < v_max){
               z_div_primme( &tmpz1, &beta, &alpha);
               H[v_size-1+v_size*v_max].r = - rnormprev/rnorm*tmpz1.r;
               H[v_size-1+v_size*v_max].i = - rnormprev/rnorm*tmpz1.i;
               H[v_size+(v_size-1)*v_max].r    =- rnorm/rnormprev*tmpz.r;
               H[v_size+(v_size-1)*v_max].i    =- rnorm/rnormprev*tmpz.i;}
       }


      stoptol_used = tol*(*normb);
      if(ConvTestOpt == 2){
        tmpc3 = wrap_zsum_cdot(&n,x,&ONE,x,&ONE,params); 
        xnorm = sqrt(tmpc3.r);
        stoptol_cur  = MACHEPS*((*AnormEst)*xnorm+(*normb));
        if(stoptol_used < stoptol_cur)
           stoptol_used = stoptol_cur;
      }
      
      
      if( rnorm < stoptol_used) { *flag = 0; printf("Stopped when norm(r) <= %g\n",stoptol_used); break;}

      /* update old parameters */
      rnormprev = rnorm;
      betaprev.r = beta.r; betaprev.i = beta.i;
      alphaprev.r = alpha.r; alphaprev.i = alpha.i;
      rhoprev.r = rho.r; rhoprev.i = rho.i;

    } /* end of the bicg iterations loop */


    *iter=it+1;

    if( (it==maxit-1) & (rnorm > tol*(*normb)) ) *flag = 1;
    if (plvl) CdisplayInfo(tol,maxit,*flag,*iter-1,reshist[it]);


    /* compute final eigenvalues and eigenvectors*/
    if(nev >0 )
       { 
           allelems = v_max * v_max;
           BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE);
           tmpi=v_size;    
 
           ZG_eval(tmpH,tmpi,v_max,Hevals,SRT_OPT,(double) epsi, Hevecsl,tmpi, Hevecsr,tmpi,&info);


           if (info != 0){ 
              fprintf(stderr, "ERROR Could not compute final eignvalues\n");
              fprintf(stderr, " info %g \n", info);
              fprintf(stderr, "Exiting...\n");
              *flag=3;
              return; }


           /* copy to single precision coefficents */
           allelems=nev*tmpi;
           for (i=0; i< allelems; i++){
               Coefr[i].r=Hevecsr[i].r;  Coefr[i].i=Hevecsr[i].i;
               Coefl[i].r=Hevecsl[i].r;  Coefl[i].i=Hevecsl[i].i;}

           Crestart_X(VL, LDVL, Coefl, n, tmpi, nev, tmp_work, tmpsize);
           Crestart_X(VR, LDVR, Coefr, n, tmpi, nev, tmp_work, tmpsize);

    
        
           for(i=0; i<nev; i++)
              {
           
                  //again, use Ap_prev to store the residual std::vector of the eigenvalue since
                  //Ap_prev space is no longer needed.

                  CcomputeResNorm(&VR[i*LDVR],&VL[i*LDVL],&lambda,&ernorm,n,Ap_prev, 
                          &xlnorm,&xrnorm,&cangle,matvec,params);

                  evals[i]=lambda; rnorms[i]=ernorm;
                  if(plvl>=3)   
                      printf("%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                      i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);
              }
       }
       


  free(H);
  free(tmpH);
  free(Hevecsl);
  free(Hevecsr);
  free(Hevecslold);
  free(Hevecsrold);
  free(Hevals);
  free(Hevalsold);
  free(Coefl);
  free(Coefr);
  free(ZCoefl); 
  free(ZCoefr); 
  return;
}
/*****************************************************************************************************/

/* EIGBICG *************************************************************************************/
void Ceigbicg_ver2(int n, int lde, Complex_C *x, Complex_C *b, float *normb, float tol, int maxit, 
                   char SRT_OPT, float epsi, int ConvTestOpt, int *iter, float *reshist, int *flag, 
                   int plvl, Complex_C *work, void (*matvec) (void *, void *, void *), 
                   void (*mathvec)(void *, void *, void *), void *params, 
                   void (*dmatvec)(void *, void *, void *), void *dparams,
                   float *AnormEst, int nev, 
                   Complex_C *evals, float *rnorms, int v_max, Complex_Z *VR, int LDVR, Complex_Z *VL, 
                   int LDVL, int esize, Complex_Z *ework)
{
  Complex_C *r, *rs, *p, *ps, *Ap, *Atps, *resid;                     
  Complex_Z alpha, alphaprev, beta, betaprev, rho, rhoprev;           
  double    rnorm, rnormprev;                                 
  int       i, j,k,it, v_size,tmpi,allelems,info,ii,jj;
  int       zs, cs, ds, tmpsize;                              
  Complex_C tmpc,tmpc1,tmpc2,tmpc3,cof1,cof2;
  Complex_Z tmpz,tmpz1;  

  Complex_Z *H,*tmpH,*Hevals,*Hevecsl, *Hevecsr;            
  Complex_Z *Hevalsold,*Hevecslold, *Hevecsrold;     
  Complex_Z *ZCoefl, *ZCoefr;              
                 

  Complex_Z  *Ap_prev, *Atps_prev, *tmp_work;                   /* pointers in ework */
  double      tmpd,tmpd1,tmpd2,ernorm,xlnorm,xrnorm, xnorm, stoptol_used, stoptol_cur;
  Complex_Z  lambda,cangle;


  float MACHEPS=1e-7;

  /* constants */
  int ONE = 1;      /* for passing 1 to BLAS routines */
  Complex_C tpone = {+1.00000000e+00,+0.0000000e00}, 
            tzero = {+0.00000000e+00,+0.0000000e00}; /*single precision Complex numbers 1 and zero */
  Complex_Z z_one = {+1.0000000000000000e+00,+0.0000000000000000e00}, 
            z_zero= {+0.0000000000000000e+00,+0.0000000000000000e00}; /*double precision Complex numbers 1 and zero */
  char cR = 'R'; char cL = 'L'; char cN ='N';                       /* Character constants */
  char cV = 'V'; char cU = 'U'; char cC ='C';

  //FILE *OFile=fopen("Ceigbicg_debug.dat","w");
  Complex_C tmpcr,tmpcl;
  float tmpfr,tmpfl;

  /* allocations */
  zs = sizeof(Complex_Z); 
  cs = sizeof(Complex_C); 
  ds = sizeof(double);


  if ((H = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
       fprintf(stderr, "ERROR Could not allocate H\n");
  if ((tmpH = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
       fprintf(stderr, "ERROR Could not allocate tmpH\n");
  if ((Hevecsl = (Complex_Z *)calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsl\n");
  if ((Hevecsr = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsr\n");
  if ((Hevecslold = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecslold\n");
  if ((Hevecsrold = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsrold\n");
  if ((Hevals = (Complex_Z *) calloc(v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevals\n");
  if ((Hevalsold = (Complex_Z *) calloc(v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevalsold\n");
  if ((ZCoefl = (Complex_Z *) calloc(2*nev*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate ZCoefl\n");
  if ((ZCoefr = (Complex_Z *) calloc(2*nev*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate ZCoefr\n");


  // Initialize H and tmpH to zeros
  for(i=0; i<v_max*v_max; i++)
     { H[i]=z_zero; tmpH[i]=z_zero;} 
 

  /* setup pointers into work */
  r    = work;
  rs   = work + lde;
  p    = work + 2*lde;
  ps   = work + 3*lde;
  Ap   = work + 4*lde;
  Atps = work + 5*lde;


  /* setup pointers into ework */
  Ap_prev  = ework;
  Atps_prev= ework+lde;
  tmp_work = ework+2*lde;
  tmpsize  = esize-2*lde;  /*leftover in tmp_work*/
  

  /*Initialization*/
  tmpc = wrap_zsum_cdot(&n, b, &ONE, b, &ONE, params); /* Norm of rhs, b */
  (*normb) = sqrt(tmpc.r);


  /* If right hand side is zero return zero solution. ITER=0 */
  if ((*normb) == 0.0) {		
     for (i=0; i<n; i++) 
             x[i]=tzero; 
     (*flag) = 0;
     (*iter) = 0; 		
     reshist[0] = 0.0;
     if (plvl) CdisplayInfo(tol,maxit,*flag,*iter,reshist[0]);
     return;
  }


  /* Zero-th residual: r = b - A*x  */
  matvec(x, r, params);			
  for (i = 0; i < n; i ++) {
      r[i].r = b[i].r - r[i].r;
      r[i].i = b[i].i - r[i].i;
   }

  /* Zero-th rs, p, and ps = Zero-th r */
  BLAS_CCOPY(&n, r, &ONE, rs, &ONE);
  BLAS_CCOPY(&n, r, &ONE, p, &ONE);
  BLAS_CCOPY(&n, rs, &ONE, ps, &ONE);

  tmpc = wrap_zsum_cdot(&n, r, &ONE, r, &ONE, params); /* Sum in double */
  reshist[0] = sqrt(tmpc.r); rnormprev=reshist[0];
  if(plvl) printf("Initial residual %g\n",reshist[0]); 

  tmpc = wrap_zsum_cdot(&n, rs, &ONE, r, &ONE, params); /* Sum in double */
  rhoprev.r =tmpc.r;  rhoprev.i=tmpc.i;
  rho=z_zero; alpha=z_zero; alphaprev=z_zero; beta=z_zero; betaprev=z_zero;



   
  v_size =0;   /* counts vectors in the eignstd::vector search spaces */
  it=0;        /* iteration counter*/

  while( it < maxit-1 ){      /*main bicg loop*/


      it = it + 1;
      matvec(p,Ap,params);
      mathvec(ps,Atps,params);

      if(v_size == v_max-1){
        for(i=0; i<n; i++)
           {
               Ap_prev[i].r   = Ap[i].r; 
               Ap_prev[i].i   = Ap[i].i;
               Atps_prev[i].r = Atps[i].r;
               Atps_prev[i].i = Atps[i].i; 
           }
       }
      
      if( (v_size == v_max) && nev >0 ){  /*eigenvalue part*/
         
         allelems=v_max*v_max;
         BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE); //copy H into a temporary storage to be used in ZG_eval  

         ZG_eval(tmpH, v_max, v_max, Hevals, SRT_OPT, (double) epsi, 
                 Hevecsl, v_max, Hevecsr, v_max, &info);


         if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of H(v_max,v_max)\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            (*flag) = 3;
            return; }


        //Estimated norm(A)
        tmpd = z_abs_primme(Hevals[v_max-1]);
        if(tmpd > (*AnormEst) )
          (*AnormEst) = tmpd;
        

         BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE);
         tmpi=v_max-1;  

         ZG_eval(tmpH,tmpi,v_max,Hevalsold,SRT_OPT, (double) epsi, Hevecslold,v_max,
                 Hevecsrold,v_max,&info);



         if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of H(v_max-1,v_max-1)\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            *flag = 3;
            return; }
          

         /* set the v_max elements of Hevecslold and Hevecsrold to zero */
         for(i=1; i <= tmpi ; i++){
           Hevecslold[i*v_max-1].r = 0.0 ; Hevecslold[i*v_max-1].i = 0.0;
           Hevecsrold[i*v_max-1].r = 0.0 ; Hevecsrold[i*v_max-1].i = 0.0;}

        /* append the nev eigenvectors computed at step (v_max-1) after the nev
           eigenvectors computed at step v_max */


        tmpi=nev*v_max;
        BLAS_ZCOPY(&tmpi, Hevecsrold, &ONE, &Hevecsr[tmpi], &ONE);
        BLAS_ZCOPY(&tmpi, Hevecslold, &ONE, &Hevecsl[tmpi], &ONE);



        /* bi-orthogonalize the nev old vectors againist the new nev vectors */
        biortho_local(Hevecsl, v_max, Hevecsr, v_max, v_max, nev+1, 2*nev, 3);
   


        /*Compute new matrix H(1:2nev,1:2nev) = Hevecl(1:v_max,1:2nev)'*H*Hevecr(1:v_max,1:2nev) */
        /* first compute H*Hevecsr and store it it Hevecsrold */

        tmpi=2*nev;
        BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &v_max, &z_one, H, &v_max,
                    Hevecsr, &v_max, &z_zero, Hevecsrold, &v_max);   
        /* next compute Hevecsl(1:v_max,1:2nev)' * Hevecsrold(1:v_max, 1:2nev) 
           and store the result in tmpH such that the leading dimension
           of tmpH will be 2*nev and the actual result will be in the 
           2nev*2nev block of tmpH */
        BLAS_ZGEMM( &cC, &cN, &tmpi, &tmpi, &v_max, &z_one, Hevecsl,
                    &v_max, Hevecsrold, &v_max, &z_zero, tmpH, &tmpi);

       /* solve the 2nev*2nev eigenvalue problem */
       /* since tmpH won't be needed afterwards, we can pass it to GEEV 
          store eigenvalues and eigenvectors in the old arrays Hevalsold,.etc*/

       ZG_eval(tmpH,tmpi,tmpi,Hevalsold,SRT_OPT,(double) epsi, 
               Hevecslold,tmpi,Hevecsrold,tmpi,&info);


       /* check that the eigenvalue problem didn't fail */
       if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of the 2nevx2nev problem\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            *flag=3;
            return; }


       /* Compute the coefficient matrices */
       BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &tmpi, &z_one, Hevecsr,
                   &v_max, Hevecsrold, &tmpi, &z_zero, ZCoefr, &v_max);

       
       BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &tmpi, &z_one, Hevecsl,
                   &v_max, Hevecslold, &tmpi, &z_zero, ZCoefl, &v_max);





       /* Copy double precision Coef matrices into single precision ones */
       //allelems=2*nev*v_max;
       //for (i=0; i< allelems; i++){
       //   Coefr[i].r=ZCoefr[i].r;  Coefr[i].i=ZCoefr[i].i;
       //   Coefl[i].r=ZCoefl[i].r;  Coefl[i].i=ZCoefl[i].i;}



       /* Restarting */
       v_size=2*nev;

       Zrestart_X(VL, LDVL, ZCoefl, n, v_max, v_size, tmp_work, tmpsize);
       Zrestart_X(VR, LDVR, ZCoefr, n, v_max, v_size, tmp_work, tmpsize);

       //-----------DEBUG---------------------------
       //Nromalize the new 2*nev vectors VL and VR:
       /*
       for(i=0; i<v_size; i++)
         {
             tmpc = wrap_zsum_cdot(&n,&VR[i*LDVR],&ONE,&VR[i*LDVR],&ONE,params);
             for(j=0; j<n; j++)
                {
                    VR[j+i*LDVR].r = VR[j+i*LDVR].r / sqrt(tmpc.r);
                    VR[j+i*LDVR].i = VR[j+i*LDVR].i / sqrt(tmpc.r);
                }
             tmpc = wrap_zsum_cdot(&n,&VR[i*LDVR],&ONE,&VL[i*LDVL],&ONE,params);
             c_div_primme(&tmpc1,&tpone,&tmpc);             
             BLAS_CSCAL(&n,&tmpc1,&VL[i*LDVL],&ONE);
         }
      */
      //-------------END DEBUG-----------------------
     



       /* Restart H = diag(Hevals) plus a column and a row */
       allelems = v_max*v_max;
       for (i = 0; i < allelems; i++ )  
             H[i]= z_zero; 

       for (i = 0; i < v_size; i++) {
           H[i*(v_max+1)].r = Hevalsold[i].r;
           H[i*(v_max+1)].i = Hevalsold[i].i;}
           

       /* Compute the elemnts 2nev+1 row and column of H */
       tmpz.r = beta.r; tmpz.i=-beta.i ;  
       BLAS_ZSCAL( &n, &tmpz, Atps_prev, &ONE);  /* beta'*Atps_prev */
       
       for(i=0; i<n; i++){
          Atps_prev[i].r = (double) Atps[i].r - Atps_prev[i].r;
          Atps_prev[i].i = (double) Atps[i].i - Atps_prev[i].i;}

       
       BLAS_ZSCAL( &n, &beta, Ap_prev, &ONE);

       for(i=0; i<n; i++){
          Ap_prev[i].r = (double) Ap[i].r - Ap_prev[i].r;
          Ap_prev[i].i = (double) Ap[i].i - Ap_prev[i].i;}

       z_div_primme(&tmpz1,&z_one,&rho);   // 1/rho

  
       for(i=0; i<v_size; i++){
          tmpz= wrap_zdot(&n, Atps_prev, &ONE, &VR[i*LDVR], &ONE, params); 
          H[i*v_max+v_size].r= rnorm*(tmpz.r*tmpz1.r-tmpz.i*tmpz1.i); 
          H[i*v_max+v_size].i= rnorm*(tmpz.r*tmpz1.i+tmpz.i*tmpz1.r);
          tmpz = wrap_zdot(&n, &VL[i*LDVL], &ONE, Ap_prev, &ONE, params);
          H[i+v_size*v_max].r = tmpz.r/rnorm; H[i+v_size*v_max].i= tmpz.i/rnorm;}

       if(plvl>=4)
         {

             for(i=0; i<2*nev; i++)
               {
                  
                  //temprarily use Ap_prev to store the eigenvalue residual std::vector.
                  //Ap_prev is not needed at this point till next restart when it will be
                  //recomputed.
                  ZcomputeResNorm(&VR[i*LDVR],&VL[i*LDVL],&lambda,&ernorm,n,Ap_prev, 
                                 &xlnorm,&xrnorm,&cangle,dmatvec,dparams);
                  printf("%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                          i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);

                  //fprintf(OFile,"%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                  //        i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);




               }
         } 


      } /*end of the eigenvalue part if(v_size==v_max & nev >0) */

      /* continue bicg from that point on */
      
      // add a new std::vector to the basis

           
      for(i=0; i<n; i++) {
          VL[LDVL*v_size + i].r = rs[i].r;  VL[LDVL*v_size+i].i = rs[i].i;
          VR[LDVR*v_size+i].r   = r[i].r ;  VR[LDVR*v_size+i].i = r[i].i ; } 
      
      tmpz.r=rnormprev; tmpz.i=0.0e+00;
      z_div_primme( &tmpz1, &z_one, &tmpz);
      BLAS_ZSCAL( &n, &tmpz1, &VR[LDVR*v_size], &ONE);

      z_div_primme( &tmpz1, &tmpz, &rhoprev);
      BLAS_ZSCAL( &n, &tmpz1, &VL[LDVL*v_size], &ONE);


      v_size++;

      tmpc = wrap_zsum_cdot(&n, ps, &ONE, Ap, &ONE, params);

      tmpz.r=tmpc.r; tmpz.i=tmpc.i;

      z_div_primme( &alpha, &rhoprev, &tmpz);  
      
      if(z_abs_primme(alpha) < 1.e-14 ) { *flag=2; break;}
      
      tmpc.r = -alpha.r; tmpc.i = -alpha.i;      //update r
      BLAS_CAXPY( &n, &tmpc, Ap, &ONE, r, &ONE);

      tmpc.r= -alpha.r; tmpc.i = alpha.i;
      BLAS_CAXPY( &n, &tmpc, Atps, &ONE, rs, &ONE);  //update rs

      tmpc = wrap_zsum_cdot(&n, r, &ONE, r, &ONE, params); /* Sum in double */
      rnorm = sqrt(tmpc.r);
      reshist[it] = rnorm;
      if(plvl>=2) printf("iter# %d, residual norm = %0.12g\n",it,rnorm);

      tmpc.r = alpha.r; tmpc.i = alpha.i;          //update the solution
      BLAS_CAXPY( &n, &tmpc, p, &ONE, x, &ONE);




      tmpc= wrap_zsum_cdot(&n, rs, &ONE, r, &ONE, params);
      rho.r = tmpc.r;  rho.i = tmpc.i;
      if(z_abs_primme(rho)< 1.e-40){ *flag=2; break;}
 
      z_div_primme( &beta, &rho, &rhoprev);
      if(z_abs_primme(beta)< 1.e-14){ *flag=2; break;}
      


      tmpc.r= beta.r; tmpc.i = beta.i;           //update p
      BLAS_CSCAL( &n, &tmpc, p, &ONE);
      BLAS_CAXPY( &n, &tpone, r, &ONE, p, &ONE);

      tmpc.r = beta.r; tmpc.i = -beta.i;         //update ps
      BLAS_CSCAL( &n, &tmpc, ps, &ONE);
      BLAS_CAXPY( &n, &tpone, rs, &ONE, ps, &ONE);


      /* elements of the projection matrix */
      /* diagonal elemnts */
      
      z_div_primme( &tmpz, &z_one, &alpha);
      if((v_size-1) == 0){ 
      H[(v_max+1)*(v_size-1)].r = tmpz.r;
      H[(v_max+1)*(v_size-1)].i = tmpz.i;}
      else{
         z_div_primme( &tmpz1, &betaprev, &alphaprev);
         H[(v_max+1)*(v_size-1)].r =  tmpz.r + tmpz1.r;
         H[(v_max+1)*(v_size-1)].i =  tmpz.i + tmpz1.i;}

      /* off-diagonal */
      if( v_size < v_max){
        z_div_primme( &tmpz1, &beta, &alpha);
        H[v_size-1+v_size*v_max].r = - rnormprev/rnorm*tmpz1.r;
        H[v_size-1+v_size*v_max].i = - rnormprev/rnorm*tmpz1.i;
        H[v_size+(v_size-1)*v_max].r    =- rnorm/rnormprev*tmpz.r;
        H[v_size+(v_size-1)*v_max].i    =- rnorm/rnormprev*tmpz.i;}


      stoptol_used = tol*(*normb);
      if(ConvTestOpt == 2){
        tmpc3 = wrap_zsum_cdot(&n,x,&ONE,x,&ONE,params); 
        xnorm = sqrt(tmpc3.r);
        stoptol_cur  = MACHEPS*((*AnormEst)*xnorm+(*normb));
        if(stoptol_used < stoptol_cur)
           stoptol_used = stoptol_cur;
      }
      
      
      if( rnorm < stoptol_used) { *flag = 0; printf("Stopped when norm(r) <= %g\n",stoptol_used); break;}

      /* update old parameters */
      rnormprev = rnorm;
      betaprev.r = beta.r; betaprev.i = beta.i;
      alphaprev.r = alpha.r; alphaprev.i = alpha.i;
      rhoprev.r = rho.r; rhoprev.i = rho.i;

    } /* end of the bicg iterations loop */


    *iter=it+1;

    if( (it==maxit-1) & (rnorm > tol*(*normb)) ) *flag = 1;


    /* compute final eigenvalues and eigenvectors*/ 
    allelems = v_max * v_max;
    BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE);
    tmpi=v_size;    
 
    ZG_eval(tmpH,tmpi,v_max,Hevals,SRT_OPT,(double) epsi, Hevecsl,tmpi, Hevecsr,tmpi,&info);


    if (info != 0){ 
      fprintf(stderr, "ERROR Could not compute final eignvalues\n");
      fprintf(stderr, " info %g \n", info);
      fprintf(stderr, "Exiting...\n");
      *flag=3;
      return; }


    /* copy to single precision coefficents */
    allelems=nev*tmpi;
    for (i=0; i< allelems; i++){
        ZCoefr[i].r=Hevecsr[i].r;  ZCoefr[i].i=Hevecsr[i].i;
        ZCoefl[i].r=Hevecsl[i].r;  ZCoefl[i].i=Hevecsl[i].i;}

    Zrestart_X(VL, LDVL, ZCoefl, n, tmpi, nev, tmp_work, tmpsize);
    Zrestart_X(VR, LDVR, ZCoefr, n, tmpi, nev, tmp_work, tmpsize);

    
    if (plvl) CdisplayInfo(tol,maxit,*flag,*iter-1,reshist[it]);
    
    for(i=0; i<nev; i++)
        {
           
           //again, use Ap_prev to store the residual std::vector of the eigenvalue since
           //Ap_prev space is no longer needed.

           ZcomputeResNorm(&VR[i*LDVR],&VL[i*LDVL],&lambda,&ernorm,n,Ap_prev, 
                          &xlnorm,&xrnorm,&cangle,dmatvec,dparams);

           evals[i].r=lambda.r;
           evals[i].i=lambda.i;
           rnorms[i]=ernorm;
           if(plvl>=3)   
               printf("%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                      i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);
        }
       


  free(H);
  free(tmpH);
  free(Hevecsl);
  free(Hevecsr);
  free(Hevecslold);
  free(Hevecsrold);
  free(Hevals);
  free(Hevalsold);
  free(ZCoefl); 
  free(ZCoefr); 
  return;
}








/**************************************************************************************************/
void Zeigbicg(int n, int lde, Complex_Z *x, Complex_Z *b, double *normb, double tol, int maxit, 
              char SRT_OPT, double epsi, int ConvTestOpt, int *iter, double *reshist, int *flag, 
              int plvl, Complex_Z *work, void (*matvec) (void *, void *, void *), 
              void (*mathvec)(void *, void *, void *), void *params, double *AnormEst, int nev, 
              Complex_Z *evals, double *rnorms, int v_max, Complex_Z *VR, int LDVR, Complex_Z *VL, 
              int LDVL, int esize, Complex_Z *ework)
{

  
  Complex_Z *r, *rs, *p, *ps, *Ap, *Atps;                     
  Complex_Z alpha, alphaprev, beta, betaprev, rho, rhoprev;           
  double    rnorm, rnormprev;                                 
  int       i, j,k,it, v_size,tmpi,allelems,info,ii,jj;
  int       zs,ds,tmpsize;                              
  Complex_Z tmpc,tmpc1,tmpc2,tmpc3,cof1,cof2;
  Complex_Z tmpz,tmpz1;  
  Complex_Z *H,*tmpH,*Hevals,*Hevecsl, *Hevecsr;            
  Complex_Z *Hevalsold,*Hevecslold, *Hevecsrold;     
  Complex_Z *ZCoefl, *ZCoefr;                               
  Complex_Z *Ap_prev,*Atps_prev,*tmp_work;                  
  double   tmpd,tmpd1,tmpd2,ernorm,xlnorm,xrnorm, xnorm, stoptol_used, stoptol_cur;
  Complex_Z lambda,cangle;

  double MACHEPS=1e-16;

  /* constants */
  int ONE = 1;      /* for passing 1 to BLAS routines */
  Complex_Z z_one = {+1.0000000000000000e+00,+0.0000000000000000e00}, 
            z_zero= {+0.0000000000000000e+00,+0.0000000000000000e00}; /*double precision Complex numbers 1 and zero */
  char cR = 'R'; char cL = 'L'; char cN ='N';                       /* Character constants */
  char cV = 'V'; char cU = 'U'; char cC ='C';

  //FILE *OFile=fopen("Zeigbicg_debug.dat","w");

  /* allocations */
  zs = sizeof(Complex_Z); 
  ds = sizeof(double);

  if ((H = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
       fprintf(stderr, "ERROR Could not allocate H\n");
  if ((tmpH = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
       fprintf(stderr, "ERROR Could not allocate tmpH\n");
  if ((Hevecsl = (Complex_Z *)calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsl\n");
  if ((Hevecsr = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsr\n");
  if ((Hevecslold = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecslold\n");
  if ((Hevecsrold = (Complex_Z *) calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsrold\n");
  if ((Hevals = (Complex_Z *) calloc(v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevals\n");
  if ((Hevalsold = (Complex_Z *) calloc(v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevalsold\n");
  if ((ZCoefl = (Complex_Z *) calloc(2*nev*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate ZCoefl\n");
  if ((ZCoefr = (Complex_Z *) calloc(2*nev*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate ZCoefr\n");


  

  // Initialize H and tmpH to zeros
  for(i=0; i<v_max*v_max; i++)
     { H[i]=z_zero; tmpH[i]=z_zero;} 
 

  /* setup pointers into work */
  r    = work;
  rs   = work + lde;
  p    = work + 2*lde;
  ps   = work + 3*lde;
  Ap   = work + 4*lde;
  Atps = work + 5*lde;


  /* setup pointers into ework */
  Ap_prev  = ework;
  Atps_prev= ework+lde;
  tmp_work = ework+2*lde;
  tmpsize  = esize-2*lde;  /*leftover in tmp_work*/
  

  /*Initialization*/
  tmpc = wrap_zdot(&n, b, &ONE, b, &ONE, params); /* Norm of rhs, b */
  (*normb) = sqrt(tmpc.r);


  /* If right hand side is zero return zero solution. ITER=0 */
  if ((*normb) == 0.0) {		
     for (i=0; i<n; i++) 
             x[i]=z_zero; 
     (*flag) = 0;
     (*iter) = 0; 		
     reshist[0] = 0.0;
     if (plvl) ZdisplayInfo(tol,maxit,*flag,*iter,reshist[0]);
     return;
  }


  /* Zero-th residual: r = b - A*x  */
  matvec(x, r, params);			
  for (i = 0; i < n; i ++) {
      r[i].r = b[i].r - r[i].r;
      r[i].i = b[i].i - r[i].i;
   }

  /* Zero-th rs, p, and ps = Zero-th r */
  BLAS_ZCOPY(&n, r, &ONE, rs, &ONE);
  BLAS_ZCOPY(&n, r, &ONE, p, &ONE);
  BLAS_ZCOPY(&n, rs, &ONE, ps, &ONE);

  tmpc = wrap_zdot(&n, r, &ONE, r, &ONE, params); /* Sum in double */
  reshist[0] = sqrt(tmpc.r); rnormprev=reshist[0];
  if(plvl) printf("Initial residual %g\n",reshist[0]); 

  tmpc = wrap_zdot(&n, rs, &ONE, r, &ONE, params); /* Sum in double */
  rhoprev.r =tmpc.r;  rhoprev.i=tmpc.i;
  rho=z_zero; alpha=z_zero; alphaprev=z_zero; beta=z_zero; betaprev=z_zero;



   
  v_size =0;   /* counts vectors in the eignstd::vector search spaces */
  it=0;        /* iteration counter*/

  while( it < maxit-1 ){      /*main bicg loop*/


      it = it + 1;

      matvec(p,Ap,params);
      mathvec(ps,Atps,params);



      if(v_size == v_max-1){
        BLAS_ZCOPY(&n, Ap, &ONE, Ap_prev, &ONE);
        BLAS_ZCOPY(&n, Atps, &ONE, Atps_prev, &ONE);}
      
      if( (v_size == v_max) && (nev >0) ){  /*eigenvalue part*/
         
         allelems=v_max*v_max;
         BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE); //copy H into a temporary storage to be used in ZG_eval  

         ZG_eval(tmpH, v_max, v_max, Hevals, SRT_OPT, epsi, 
                 Hevecsl, v_max, Hevecsr, v_max, &info);


         if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of H(v_max,v_max)\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            (*flag) = 3;
            return; }


        //Estimated norm(A)
        tmpd = z_abs_primme(Hevals[v_max-1]);
        if(tmpd > (*AnormEst) )
          (*AnormEst) = tmpd;
        

         BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE);
         tmpi=v_max-1;  

         ZG_eval(tmpH,tmpi,v_max,Hevalsold,SRT_OPT, epsi, Hevecslold,v_max,
                 Hevecsrold,v_max,&info);



         if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of H(v_max-1,v_max-1)\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            *flag = 3;
            return; }
          

         /* set the v_max elements of Hevecslold and Hevecsrold to zero */
         for(i=1; i <= tmpi ; i++){
           Hevecslold[i*v_max-1].r = 0.0 ; Hevecslold[i*v_max-1].i = 0.0;
           Hevecsrold[i*v_max-1].r = 0.0 ; Hevecsrold[i*v_max-1].i = 0.0;}

        /* append the nev eigenvectors computed at step (v_max-1) after the nev
           eigenvectors computed at step v_max */


        tmpi=nev*v_max;
        BLAS_ZCOPY(&tmpi, Hevecsrold, &ONE, &Hevecsr[tmpi], &ONE);
        BLAS_ZCOPY(&tmpi, Hevecslold, &ONE, &Hevecsl[tmpi], &ONE);



        /* bi-orthogonalize the nev old vectors againist the new nev vectors */
        biortho_local(Hevecsl, v_max, Hevecsr, v_max, v_max, nev+1, 2*nev, 3);
   


        /*Compute new matrix H(1:2nev,1:2nev) = Hevecl(1:v_max,1:2nev)'*H*Hevecr(1:v_max,1:2nev) */
        /* first compute H*Hevecsr and store it it Hevecsrold */

        tmpi=2*nev;
        BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &v_max, &z_one, H, &v_max,
                    Hevecsr, &v_max, &z_zero, Hevecsrold, &v_max);   
        /* next compute Hevecsl(1:v_max,1:2nev)' * Hevecsrold(1:v_max, 1:2nev) 
           and store the result in tmpH such that the leading dimension
           of tmpH will be 2*nev and the actual result will be in the 
           2nev*2nev block of tmpH */
        BLAS_ZGEMM( &cC, &cN, &tmpi, &tmpi, &v_max, &z_one, Hevecsl,
                    &v_max, Hevecsrold, &v_max, &z_zero, tmpH, &tmpi);

       /* solve the 2nev*2nev eigenvalue problem */
       /* since tmpH won't be needed afterwards, we can pass it to GEEV 
          store eigenvalues and eigenvectors in the old arrays Hevalsold,.etc*/

       ZG_eval(tmpH,tmpi,tmpi,Hevalsold,SRT_OPT,epsi, 
               Hevecslold,tmpi,Hevecsrold,tmpi,&info);


       /* check that the eigenvalue problem didn't fail */
       if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute eignvalues of the 2nevx2nev problem\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            *flag=3;
            return; }


       /* Compute the coefficient matrices */
       BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &tmpi, &z_one, Hevecsr,
                   &v_max, Hevecsrold, &tmpi, &z_zero, ZCoefr, &v_max);

       
       BLAS_ZGEMM( &cN, &cN, &v_max, &tmpi, &tmpi, &z_one, Hevecsl,
                   &v_max, Hevecslold, &tmpi, &z_zero, ZCoefl, &v_max);


       /* Restarting */
       v_size=2*nev;

       Zrestart_X(VL, LDVL, ZCoefl, n, v_max, v_size, tmp_work, tmpsize);
       Zrestart_X(VR, LDVR, ZCoefr, n, v_max, v_size, tmp_work, tmpsize);

       /* Restart H = diag(Hevals) plus a column and a row */
       allelems = v_max*v_max;
       for (i = 0; i < allelems; i++ )  
             H[i]= z_zero; 

       for (i = 0; i < v_size; i++) {
           H[i*(v_max+1)].r = Hevalsold[i].r;
           H[i*(v_max+1)].i = Hevalsold[i].i;}
           

       /* Compute the elemnts 2nev+1 row and column of H */
       tmpc.r = beta.r; tmpc.i=-beta.i ;  /* beta' to a single precision: WARNING */
       BLAS_ZSCAL( &n, &tmpc, Atps_prev, &ONE);  /* beta'*Atps_prev */
       
       for(i=0; i<n; i++){
          Atps_prev[i].r = Atps[i].r - Atps_prev[i].r;
          Atps_prev[i].i = Atps[i].i - Atps_prev[i].i;}

       tmpc.r = beta.r; tmpc.i= beta.i;
       BLAS_ZSCAL( &n, &tmpc, Ap_prev, &ONE);

       for(i=0; i<n; i++){
          Ap_prev[i].r = Ap[i].r - Ap_prev[i].r;
          Ap_prev[i].i = Ap[i].i - Ap_prev[i].i;}

       z_div_primme(&tmpz1,&z_one,&rho);   // 1/rho

  
       for(i=0; i<v_size; i++){
          tmpc= wrap_zdot(&n, Atps_prev, &ONE, &VR[i*LDVR], &ONE, params); 
          H[i*v_max+v_size].r= rnorm*(tmpc.r*tmpz1.r-tmpc.i*tmpz1.i); 
          H[i*v_max+v_size].i= rnorm*(tmpc.r*tmpz1.i+tmpc.i*tmpz1.r);
          tmpc = wrap_zdot(&n, &VL[i*LDVL], &ONE, Ap_prev, &ONE, params);
          H[i+v_size*v_max].r = tmpc.r/rnorm; H[i+v_size*v_max].i= tmpc.i/rnorm;}

       if(plvl>=4)
         {

             for(i=0; i<2*nev; i++)
               {
                  
                  //temprarily use Ap_prev to store the eigenvalue residual std::vector.
                  //Ap_prev is not needed at this point till next restart when it will be
                  //recomputed.
                  ZcomputeResNorm(&VR[i*LDVR],&VL[i*LDVL],&lambda,&ernorm,n,Ap_prev, 
                                 &xlnorm,&xrnorm,&cangle,matvec,params);
                  printf("%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                          i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);

                  //fprintf(OFile,"%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                  //        i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);




               }
         } 


      } /*end of the eigenvalue part if(v_size==v_max & nev >0) */

      /* continue bicg from that point on */
      
      // add a new std::vector to the basis

      if(nev>0)
        {     
            for(i=0; i<n; i++) {
                VL[LDVL*v_size + i].r = rs[i].r;  VL[LDVL*v_size+i].i = rs[i].i;
                VR[LDVR*v_size+i].r   = r[i].r ;  VR[LDVR*v_size+i].i = r[i].i ; } 
      
             tmpz.r=rnormprev; tmpz.i=0.0e+00;
             z_div_primme( &tmpz1, &z_one, &tmpz);
             tmpc.r = tmpz1.r;   tmpc.i = tmpz1.i;
             BLAS_ZSCAL( &n, &tmpc, &VR[LDVR*v_size], &ONE);

             z_div_primme( &tmpz1, &tmpz, &rhoprev);
             tmpc.r= tmpz1.r;  tmpc.i = - tmpz1.i;
             BLAS_ZSCAL( &n, &tmpc, &VL[LDVL*v_size], &ONE);

            v_size++;
      }

 
      tmpc = wrap_zdot(&n, ps, &ONE, Ap, &ONE, params);

      tmpz.r=tmpc.r; tmpz.i=tmpc.i;


      z_div_primme( &alpha, &rhoprev, &tmpz);  
     
 
      if(z_abs_primme(alpha) < 1.e-14 ) { *flag=2; break;}
      
      tmpc.r = -alpha.r; tmpc.i = -alpha.i;      //update r
      BLAS_ZAXPY( &n, &tmpc, Ap, &ONE, r, &ONE);

      tmpc.r= -alpha.r; tmpc.i = alpha.i;
      BLAS_ZAXPY( &n, &tmpc, Atps, &ONE, rs, &ONE);  //update rs

      tmpc = wrap_zdot(&n, r, &ONE, r, &ONE, params); /* Sum in double */
      rnorm = sqrt(tmpc.r);
      reshist[it] = rnorm;
      if(plvl>=2) printf("iter# %d, residual norm = %0.12g\n",it,rnorm);

      tmpc.r = alpha.r; tmpc.i = alpha.i;          //update the solution
      BLAS_ZAXPY( &n, &tmpc, p, &ONE, x, &ONE);

      tmpc= wrap_zdot(&n, rs, &ONE, r, &ONE, params);
      rho.r = tmpc.r;  rho.i = tmpc.i;
      if(z_abs_primme(rho)< 1.e-40){ *flag=2; break;}
 
      z_div_primme( &beta, &rho, &rhoprev);
      if(z_abs_primme(beta)< 1.e-14){ *flag=2; break;}
      


      tmpc.r= beta.r; tmpc.i = beta.i;           //update p
      BLAS_ZSCAL( &n, &tmpc, p, &ONE);
      BLAS_ZAXPY( &n, &z_one, r, &ONE, p, &ONE);

      tmpc.r = beta.r; tmpc.i = -beta.i;         //update ps
      BLAS_ZSCAL( &n, &tmpc, ps, &ONE);
      BLAS_ZAXPY( &n, &z_one, rs, &ONE, ps, &ONE);


      /* elements of the projection matrix */
      /* diagonal elemnts */
      
     
      if(nev>0)
        {
            z_div_primme( &tmpz, &z_one, &alpha);
            if((v_size-1) == 0){ 
               H[(v_max+1)*(v_size-1)].r = tmpz.r;
               H[(v_max+1)*(v_size-1)].i = tmpz.i;}
            else{
               z_div_primme( &tmpz1, &betaprev, &alphaprev);
               H[(v_max+1)*(v_size-1)].r =  tmpz.r + tmpz1.r;
               H[(v_max+1)*(v_size-1)].i =  tmpz.i + tmpz1.i;}

            /* off-diagonal */
            if( v_size < v_max){
                   z_div_primme( &tmpz1, &beta, &alpha);
                  H[v_size-1+v_size*v_max].r = - rnormprev/rnorm*tmpz1.r;
                  H[v_size-1+v_size*v_max].i = - rnormprev/rnorm*tmpz1.i;
                  H[v_size+(v_size-1)*v_max].r    =- rnorm/rnormprev*tmpz.r;
                  H[v_size+(v_size-1)*v_max].i    =- rnorm/rnormprev*tmpz.i;}
        }

      stoptol_used = tol*(*normb);
      if(ConvTestOpt == 2){
        tmpc3 = wrap_zdot(&n,x,&ONE,x,&ONE,params); 
        xnorm = sqrt(tmpc3.r);
        stoptol_cur  = MACHEPS*((*AnormEst)*xnorm+(*normb));
        if(stoptol_used < stoptol_cur)
           stoptol_used = stoptol_cur;
      }
      
      
      if( rnorm < stoptol_used) { *flag = 0; printf("Stopped when norm(r) <= %g\n",stoptol_used); break;}

      /* update old parameters */
      rnormprev = rnorm;
      betaprev.r = beta.r; betaprev.i = beta.i;
      alphaprev.r = alpha.r; alphaprev.i = alpha.i;
      rhoprev.r = rho.r; rhoprev.i = rho.i;

    } /* end of the bicg iterations loop */


    *iter=it+1;

    if( (it==maxit-1) & (rnorm > tol*(*normb)) ) *flag = 1;

    if (plvl) ZdisplayInfo(tol,maxit,*flag,*iter-1,reshist[it]);

    /* compute final eigenvalues and eigenvectors*/
    if(nev>0)
      { 
         allelems = v_max * v_max;
         BLAS_ZCOPY(&allelems, H, &ONE, tmpH, &ONE);
         tmpi=v_size;    
 
         ZG_eval(tmpH,tmpi,v_max,Hevals,SRT_OPT,epsi, Hevecsl,tmpi, Hevecsr,tmpi,&info);


         if (info != 0){ 
            fprintf(stderr, "ERROR Could not compute final eignvalues\n");
            fprintf(stderr, " info %g \n", info);
            fprintf(stderr, "Exiting...\n");
            *flag=3;
         return; }

         Zrestart_X(VL, LDVL, Hevecsl, n, tmpi, nev, tmp_work, tmpsize);
         Zrestart_X(VR, LDVR, Hevecsr, n, tmpi, nev, tmp_work, tmpsize);

    
    
         for(i=0; i<nev; i++)
            {
           
                //again, use Ap_prev to store the residual std::vector of the eigenvalue since
                //Ap_prev space is no longer needed.

                ZcomputeResNorm(&VR[i*LDVR],&VL[i*LDVL],&lambda,&ernorm,n,Ap_prev, 
                          &xlnorm,&xrnorm,&cangle,matvec,params);

                evals[i]=lambda; rnorms[i]=ernorm;
                if(plvl>=3)   
                   printf("%d, eval %g %g, ernorm %g, xlnorm %g, xrnorm %g, cangle %g %g\n",
                      i,lambda.r,lambda.i,ernorm,xlnorm,xrnorm,cangle.r,cangle.i);
           }
       
        }

  free(H);
  free(tmpH);
  free(Hevecsl);
  free(Hevecsr);
  free(Hevecslold);
  free(Hevecsrold);
  free(Hevals);
  free(Hevalsold);
  free(ZCoefl); 
  free(ZCoefr); 
  return;
}
/*****************************************************************************************************/








/********************************************************************************************/
void CcomputeResNorm( Complex_C *xr, Complex_C *xl, Complex_C *lambda, float *rnorm, int n, Complex_C *Res, 
      float *xlnorm, float *xrnorm, Complex_C *cangle, void (*matvec)(void *, void *, void *), void *params)
{
   int i, ONE = 1;
   Complex_C tempc,tempc1,tpone={1.000e+00,0.0000e+00};
   float normx2;
   
   /* Norm of left eigenstd::vector */
   tempc = wrap_zsum_cdot(&n, xl, &ONE, xl, &ONE, params);
   (*xlnorm) = sqrt(tempc.r);

   /* Norm of x squared */
   tempc1 = wrap_zsum_cdot(&n, xr, &ONE, xr, &ONE, params); /* Sum in double */
   (*xrnorm) = sqrt(tempc1.r);
   normx2 = tempc1.r; 

   /* cos of the angle between left and right eigenstd::vector given by cangle= (xl'*xr) / (norm(xl)*norm(xr)) */
   tempc = wrap_zsum_cdot(&n, xl, &ONE, xr, &ONE, params);
   (*cangle).r = tempc.r/(*xlnorm)/(*xrnorm);
   (*cangle).i = tempc.i/(*xlnorm)/(*xrnorm);


   /* compute Ax */
   matvec(xr,Res,params);

   /* lambda = xl'Axr/xl'xr */
   tempc1 = wrap_zsum_cdot(&n, xl, &ONE, Res, &ONE, params); /* Sum in double */
   c_div_primme(lambda,&tempc1,&tempc);

   /* res = w - lambda v */
   for (i=0; i<n; i++) {
        c_mul_primme(&tempc,lambda,&xr[i]);
        Res[i].r = Res[i].r-tempc.r;
        Res[i].i = Res[i].i-tempc.i;
   }
   tempc = wrap_zsum_cdot(&n, Res, &ONE, Res, &ONE, params); /* Sum in single */
   (*rnorm) = sqrt(tempc.r)/(*xrnorm);
}

/********************************************************************************************/


/******************************************************************************************/
void ZcomputeResNorm( Complex_Z *xr, Complex_Z *xl, Complex_Z *lambda, double *rnorm, int n, Complex_Z *Res, 
      double *xlnorm, double *xrnorm, Complex_Z *cangle, void (*matvec)(void *, void *, void *), void *params)
{
   int i, ONE = 1;
   Complex_Z tempc,tempc1,tpone={1.000e+00,0.0000e+00};
   double normx2;
   
   /* Norm of left eigenstd::vector */
   tempc = wrap_zdot(&n, xl, &ONE, xl, &ONE, params);
   (*xlnorm) = sqrt(tempc.r);

   /* Norm of x squared */
   tempc1 = wrap_zdot(&n, xr, &ONE, xr, &ONE, params); /* Sum in double */
   (*xrnorm) = sqrt(tempc1.r);
   normx2 = tempc1.r; 

   /* cos of the angle between left and right eigenstd::vector given by cangle= (xl'*xr) / (norm(xl)*norm(xr)) */
   tempc = wrap_zdot(&n, xl, &ONE, xr, &ONE, params);
   (*cangle).r = tempc.r/(*xlnorm)/(*xrnorm);
   (*cangle).i = tempc.i/(*xlnorm)/(*xrnorm);


   /* compute Ax */
   matvec(xr,Res,params);

   /* lambda = xl'Axr/xl'xr */
   tempc1 = wrap_zdot(&n, xl, &ONE, Res, &ONE, params); /* Sum in double */
   z_div_primme(lambda,&tempc1,&tempc);

   /* res = w - lambda v */
   for (i=0; i<n; i++) {
        z_mul_primme(&tempc,lambda,&xr[i]);
        Res[i].r = Res[i].r-tempc.r;
        Res[i].i = Res[i].i-tempc.i;
   }
   tempc = wrap_zdot(&n, Res, &ONE, Res, &ONE, params); /* Sum in single */
   (*rnorm) = sqrt(tempc.r)/(*xrnorm);
}
/************************************************************************************************************/

/************************************************************************************************************/
static void CdisplayInfo(float tol, int maxit, int flag, int iter, float resnorm) 
{
  if (flag != 0) {
    printf("BICG stopped at iteration %d with flag %d. ", iter, flag);
  }
  
  switch(flag) {
  case 0:
    if (iter == 0)
      printf("The initial guess has relative residual %0.5g which is within\nthe desired tolerance %0.5g\n", resnorm, tol);
    else
      printf("BICG converged at iteration %d to a solution with residual norm %0.5g", iter, resnorm);
    break;
  case 1:
    printf("\nbecause the maximum number of iterations was reached.");
    break;
  case 2:
    printf("\nbecause a scalar quantity became too small.");
    break;
  }
  
  if (flag != 0)
    printf("\nThe iterate returned at iteration %d has residual norm %0.5g",iter,resnorm);

  printf("\n");
}
/************************************************************************************************************/
static void ZdisplayInfo(double tol, int maxit, int flag, int iter, double resnorm) 
{
  if (flag != 0) {
    printf("BICG stopped at iteration %d with flag %d. ", iter, flag);
  }
  
  switch(flag) {
  case 0:
    if (iter == 0)
      printf("The initial guess has relative residual %0.5g which is within\nthe desired tolerance %0.5g\n", resnorm, tol);
    else
      printf("BICG converged at iteration %d to a solution with residual norm %0.5g", iter, resnorm);
    break;
  case 1:
    printf("\nbecause the maximum number of iterations was reached.");
    break;
  case 2:
    printf("\nbecause a scalar quantity became too small.");
    break;
  }
  
  if (flag != 0)
    printf("\nThe iterate returned at iteration %d has residual norm %0.5g",iter,resnorm);

  printf("\n");
}




