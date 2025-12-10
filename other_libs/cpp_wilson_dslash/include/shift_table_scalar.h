#ifndef SHIFT_TABLE_SCALAR_H
#define SHIFT_TABLE_SCALAR_H

#include <cstdlib>
#include <cache.h>


namespace CPlusPlusWilsonDslash {

  class ShiftTable {
  public:

    ShiftTable(const int latt_size[],
	       void (*getSiteCoords)(int coord[], int node, int linearsite),
	       int (*getLinearSiteIndex)(const int coord[]),
	       int (*nodeNum)(const int coord[])
	       );
    
    ~ShiftTable() {
      free(xshift_table);
      free(xsite_table);
      //      free(xpath_table);
    }

    inline
    int siteTable(int i) {
      return site_table[i];
    }

    inline
    int forwardNeighbor(int mysite, int mymu) {
      return shift_table[mymu + 4*(mysite + total_vol)];
    }

    inline
    int backwardNeighbor(int mysite, int mymu) {
      return shift_table[mymu + 4*mysite];
    }
    
    inline int totalVolCB() {
      return total_vol_cb;
    }

    //  /* my_index = cb*total_vol_cb + site */
    // inline int getPathSite(int my_index) const {
    //   return path_table[my_index];
    // }



  private:
    /* Tables */
    int *xshift_table;        /* Unaligned */
    int *shift_table;         /* Aligned */

    int *xsite_table;         /* Unaligned */
    int *site_table;          /* Aligned */

    //    int *xpath_table;         /* Unaligned */
    // int *path_table;          /* Aligned */
   

    int tot_size[4];          /* Class scope members */
    int total_vol;
    int total_vol_cb;         /* Useful numbers */
    const int Nd;                   /* No of Dimensions */
    
    /* Functions for making the shift tables */
    void mySiteCoords4D(int gcoords[], int node, int linearsite);
    int myLinearSiteIndex4D(const int gcoords[]);



   
    
    inline
    int localSite4d(int coord[], int latt_size[])
    {
      int order = 0;
      int mmu;
      
      for(mmu=4-1; mmu >= 1; --mmu) {
	order = latt_size[mmu-1]*(coord[mmu] + order);
      }
      order += coord[0];
      
      return order;
    }

    inline 
    void crtesn4d(int ipos, const int latt_size[], int coord[] )
    {
      int Ndim=0; /* Start running x fastest */
      int i, ix;
      
      /* Calculate the Cartesian coordinates of the VALUE of IPOS where the 
      * value is defined by
      *
      *     for i = 0 to NDIM-1  {
      *        X_i  <- mod( IPOS, L(i) )
      *        IPOS <- int( IPOS / L(i) )
      *     }
      *
      * NOTE: here the coord(i) and IPOS have their origin at 0. 
      */
      for(i = Ndim; i < Ndim+4; ++i) {
	ix=i%4;  /* This lets me start with the time direction and then wraparound */
	
	coord[ix] = ipos % latt_size[ix];
	ipos = ipos / latt_size[ix];
      }
    }

    inline 
    void offs(int temp[], const int coord[], int mu, int isign)
    {
      int i;
      
      for(i=0; i < 4; ++i) {
	temp[i] = coord[i];
      }
      
      /* translate address to neighbour */
      temp[mu] = (temp[mu] + isign + 2*tot_size[mu]) % tot_size[mu];
    } 

    inline
    int parity(const int coord[])
    {
      int m;
      int sum = 0;
      
      for(m=0; m < 4; ++m) {
	sum += coord[m];
      }
      
      return sum % 2;
    }



#if 1
    void setupPathTable(int (*getLinearSiteIndex)(const int coord[]),
			int* inv_table) ;    
#endif

  };


  
} // Namespace

#endif
