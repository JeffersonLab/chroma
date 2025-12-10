/*------------------------------------------------------------------------
  Index sorting of a ingle or double precision array 
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
--------------------------------------------------------------------------*/

#include "sort.h"


void index_bubble_sort(unsigned long n, double *arr, unsigned long *indx)
{

   unsigned long indxt;
   int i, len, isSorted = 0;
   
   /* Intialize the index array */ 
   for( i=0 ; i<n ; i++) indx[i]=i;

   /* Perform the bubble sort on the index array */
   for( len = n-1; len> 0; len--)
   {
      isSorted = 1;
      for(i=0; i< len; i++){
        if(arr[indx[i]] > arr[indx[i+1]] )
           {
             isSorted = 0; /* Not finished yet */
             /* swap the indeces */
             indxt = indx[i];
             indx[i] = indx[i+1];
             indx[i+1] = indxt ;
           }
       }
       if(isSorted == 1) break;
      
   } 
 
}




void index_bubble_sort_single(unsigned long n, float *arr, unsigned long *indx)
{

   unsigned long indxt;
   int i, len, isSorted = 0;
   
   /* Intialize the index array */ 
   for( i=0 ; i<n ; i++) indx[i]=i;

   /* Perform the bubble sort on the index array */
   for( len = n-1; len> 0; len--)
   {
      isSorted = 1;
      for(i=0; i< len; i++){
        if(arr[indx[i]] > arr[indx[i+1]] )
           {
             isSorted = 0; /* Not finished yet */
             /* swap the indeces */
             indxt = indx[i];
             indx[i] = indx[i+1];
             indx[i+1] = indxt ;
           }
       }
       if(isSorted == 1) break;
      
   } 
 
}



