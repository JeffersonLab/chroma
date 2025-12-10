/*------------------------------------------------------------------------
  Index sorting of a single or double precision array 
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
--------------------------------------------------------------------------*/
#ifndef SORT_H_
#define SORT_H_



void index_bubble_sort(unsigned long N, double *arr, unsigned long *indx);
/* Indexes an array arr[0..n-1], i.e.,outputs the array indx[0..n-1] such that 
arr[indx[j]] is in ascending order for j = 0,1,..,n-1. The input quantities n
and arr are not changed.
-------------------------------------------------------------------------------------------*/

void index_bubble_sort_single(unsigned long N, float *arr, unsigned long *indx);
/*same as index_bubble_sort for float */

#endif
