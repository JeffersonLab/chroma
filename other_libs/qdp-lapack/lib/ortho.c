/* $Id: ortho.c,v 1.2 2008-04-22 03:53:05 kostas Exp $ */
/**********************************************************************
 * Function ortho_retained_vectors -- This function orthogonalizes
 *   coefficient vectors (the eigenvectors of the projection H) that
 *   are retained from a previous iteration.  The retained coefficient
 *   vectors are orthogonalized versus the coefficient vectors of the
 *   restarted basis.  This orthogonalizes previous and current Ritz 
 *   vectors cheaply and without communication because the coefficient 
 *   vectors are stored by each process.
 *   
 * Input parameters
 * ----------------
 * OrthogonalVecs  Coefficient vectors from the current iteration to 
 *    be orthogonalized against.
 *
 * length          The dimension of each std::vector
 *
 * numVectors      The number of current vectors
 * 
 * 
 * Input/Output parameters
 * -----------------------
 * newVectors The coefficient vectors obtained from the previous iteration
 * 		   Note that their leading dimension is primme->maxBasisSize.
 *
 * numNew     The number of previous vectors to retain
 *
 * rwork           work array of size max(numVectors, numNew).  If it
 *                 is of size maxBasisSize, then it will always be sufficiently
 *                 large.
 *
 *
 * Return value
 * ------------
 * int   The number of previous vectors successfully orthogonalized
 *
 ****************************************************************************/
#include <math.h>
#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical.h"
#include "qdp-lapack_numerical_private.h"

int ortho_new_vectors (Complex_C *OrthogonalVecs, int ldV, int length, 
  int numVectors, int numReorthos, Complex_C *newVectors, int numNew, 
  double machEps, Complex_C *rwork, Complex_C *work, void *params ) {

   int i;       /* Loop counter                                     */
   int tmpsize;
   int nOrths;  /* Number of times a std::vector has been orthogonalized */
   int zeroed;  /* True if the std::vector norm was reduced below 1e-14  */
   int ONE = 1;
   double norm; /* Vector norm.                                     */
   Complex_C ztmp;/* Temp accumulation var  			    */
   		/* and some constants				    */
   Complex_C tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
   Complex_C tmone = {-1.0e+00,+0.0e00};
   char cR = 'R'; char cL = 'L'; char cN ='N';
   char cV = 'V'; char cU = 'U'; char cC ='C';


   /* Orthogonalize each of the numNew vectors against the current */
   /* vectors and amongst themselves.                                   */

   i = 0;  

   while (i < numNew) {
      zeroed = 0;

      /* Orthogonalize each std::vector twice to ensure numerical stability */
      /* Still I believe orthogonalizing only once may be ok ???????*/

      for (nOrths = 0; nOrths < numReorthos; nOrths++) {

         /* Orthogonalize versus the numVectors current vectors */

         wrap_cgemv(&cC, &length, &numVectors, &tpone, OrthogonalVecs, 
	    &ldV, &newVectors[ldV*i], &ONE, &tzero, rwork, &ONE, 
	    work, params);

         BLAS_CGEMV(&cN, &length, &numVectors, &tmone, OrthogonalVecs, 
	    &ldV, rwork, &ONE, &tpone, &newVectors[ldV*i], &ONE);

         /* Orthogonalize against the i previous vectors that have */
         /* been orthogonalized thus far.                          */

         if (i > 0) {
            wrap_cgemv(&cC, &length, &i, &tpone, newVectors, 
	       &ldV, &newVectors[ldV*i], &ONE, &tzero, rwork, &ONE,
	       work, params);

            BLAS_CGEMV(&cN, &length, &i, &tmone, newVectors, 
	       &ldV, rwork, &ONE, &tpone, &newVectors[ldV*i], &ONE);
	 }

         ztmp = wrap_cdot(&length, &newVectors[ldV*i], &ONE,
            			   &newVectors[ldV*i], &ONE, params);
	 norm = ztmp.r;
         norm = sqrt(norm);

         if (norm < 5.0L*machEps) {
            numNew--;
            tmpsize = ldV*(numNew-i);
            BLAS_CCOPY(&tmpsize,
	      &newVectors[ldV*(i+1)], &ONE, &newVectors[ldV*i], &ONE);
            zeroed = 1;
            break;
         }

	 {ztmp.r = 1.0L/norm; ztmp.i = 0.0L;}
         BLAS_CSCAL(&length, &ztmp, &newVectors[ldV*i], &ONE);

      } /* for as many reorthos */

      if (!zeroed) {
         i++;
      }

   } /* main while loop */

   return numNew;
}
