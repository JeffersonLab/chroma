#ifndef SHIFT_TABLE_3D_SCALAR_H
#define SHIFT_TABLE_3D_SCALAR_H


#include <iostream>
#include <cache.h>
#include <cstdlib>
#include "allocate.h"

namespace CPlusPlusWilsonDslash {

  class ShiftTable3D {
  public:

    
    ~ShiftTable3D() {
      CPlusPlusWilsonDslash::dealloc(xshift_table);
      CPlusPlusWilsonDslash::dealloc(xsite_table);
    }

    inline
    int siteTable(int i) {
      return site_table[i];
    }

    inline
    int forwardNeighbor(int mysite, int mymu) {
      return shift_table[mymu + 3*(mysite + total_vol)];
    }

    inline
    int backwardNeighbor(int mysite, int mymu) {
      return shift_table[mymu + 3*mysite];
    }
    
    inline int totalVolCB() {
      return total_vol_cb;
    }
  private:
    /* Tables */
    int *xshift_table;        /* Unaligned */
    int *shift_table;         /* Aligned */
    
    int *xsite_table;         /* Unaligned */
    int *site_table;          /* Aligned */
    
    
    int tot_size[4];          /* Class scope members */
    int total_vol;
    int total_vol_cb;         /* Useful numbers */

    const int Nd3;                   /* No of Dimensions */
    
    /* Functions for making the shift tables */
    void mySiteCoords3D(int gcoords[], int node, int linearsite) 
    {
      int mu;
      int subgrid_cb_nrow[4];
      int tmp_coord[4];
      int cb3,cbb3;

      for(mu=0; mu < 4; mu++) { 
	subgrid_cb_nrow[mu] = tot_size[mu];
      }
      subgrid_cb_nrow[0] /=2;  /* Checkerboarding */
      
      /* Single processor -- all coords 0 */
      for(mu=0; mu < 4; mu++) { 
	gcoords[mu] = 0;
      }
    
      cb3=linearsite/total_vol_cb;
      
      crtesn(linearsite % total_vol_cb, subgrid_cb_nrow, tmp_coord);
      
      // Add on position within the node
      // NOTE: the cb for the x-coord is not yet determined
      gcoords[0] += 2*tmp_coord[0];
      for(mu=1; mu < 4; ++mu) {
	gcoords[mu] += tmp_coord[mu];
      }
      
      cbb3 = cb3;
      for(mu=1; mu < 3; ++mu) {
	cbb3 += gcoords[mu];
      }
      gcoords[0] += (cbb3 & 1);
    }
    
    int myLinearSiteIndex3D(const int gcoords[])
    {
      int mu;
      int subgrid_cb_nrow[4];
      int subgrid_cb_coord[4];
      int cb3;

      for(mu=0; mu < 4; mu++) { 
	subgrid_cb_nrow[mu] = tot_size[mu];
      }
      subgrid_cb_nrow[0] /=2;  /* Checkerboarding */

      cb3=0;
      for(mu=0; mu < 3; ++mu) { 
	cb3 += gcoords[mu];
      }
      cb3 &=1;
    
      subgrid_cb_coord[0] = (gcoords[0]/2)% subgrid_cb_nrow[0];
      for(mu=1; mu < 4; mu++) { 
	subgrid_cb_coord[mu] = gcoords[mu] % subgrid_cb_nrow[mu];
      }

      return localSite(subgrid_cb_coord, subgrid_cb_nrow) + cb3*total_vol_cb;
    }

    
    inline
    int localSite(int coord[], int latt_size[])
    {  
      int order = 0;
      int mmu;
      
      // In the 4D Case: t+Lt(x + Lx(y + Ly*z)
      // essentially  starting from i = dim[Nd-2]
      //  order =  latt_size[i-1]*(coord[i])
      //   and need to wrap i-1 around to Nd-1 when it gets below 0
      for(mmu=2; mmu >= 0; --mmu) {
	int wrapmu = (mmu-1) % 4;
	if ( wrapmu < 0 ) wrapmu += 4;
	order = latt_size[wrapmu]*(coord[mmu] + order);
      }
      
      order += coord[ 3 ]; /* T is fastest running */
      
      return order;
    }

    inline 
    void crtesn(int ipos, const int latt_size[], int coord[] )
    {   
      int i, ix;
      int Nd3=3;

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
      for(i = Nd3; i < Nd3+4; ++i) {
	ix=i%4;  /* This lets me start with the time direction and then wraparound */
	
	coord[ix] = ipos % latt_size[ix];
	ipos = ipos / latt_size[ix];
      }
    }

    inline 
    void offs(int temp[], const int coord[], int mu, int isign)
    {
      int i;
      
      // All for dims
      for(i=0; i < 4; ++i)
	temp[i] = coord[i];
      
      /* translate address to neighbour */
      temp[mu] = (temp[mu] + isign + 2*tot_size[mu]) % tot_size[mu];
    }
   

    inline
    int parity(const int coord[])
    {
      int m;
      int sum = 0;
      
      for(m=0; m < 3; ++m) {
	sum += coord[m];
      }
      
      return sum % 2;
    }

  public:
    ShiftTable3D(const int latt_size[],
			       void (*getSiteCoords)(int coord[], int node, int linearsite),
			       int (*getLinearSiteIndex)(const int coord[]),
			       int (*nodeNum)(const int coord[])
			       ) : Nd3(3)
    { 
      int dir; 

      int linear;
      int backward=0;
      int forward =1;
      int cb3;
      

      /* Set the lattice size, get total volume and checkerboarded volume */
      total_vol = 1;
      
      for(int mu=0; mu < 4; mu++) {
	tot_size[mu] = latt_size[mu];
	total_vol *= latt_size[mu];
      }
      total_vol_cb = total_vol / 2;

      /* Allocate the shift table */
      /* Nd3 directions
	 total vol sites
	 2 types (FWD, BKWD)
      */
      xshift_table = (int *)CPlusPlusWilsonDslash::alloc(3*total_vol*2*sizeof(int) + Cache::CacheLineSize);
      if ( xshift_table == 0x0 ) { 
	std::cerr << "Couldnt allocate xshift table " << std::endl;
	exit(1);
      }
      
      ptrdiff_t pad = 0;
      if( (ptrdiff_t)xshift_table % Cache::CacheLineSize != 0 ) { 
	pad = Cache::CacheLineSize - (ptrdiff_t)xshift_table % Cache::CacheLineSize;
      }
      shift_table = (int *)((unsigned char *)xshift_table + pad);

      
      /* Allocate the site table and the shift table */
      /* Now I want to build the site table */
      /* I want it cache line aligned? */
      xsite_table = (int *)CPlusPlusWilsonDslash::alloc(sizeof(int)*total_vol+Cache::CacheLineSize);
      if(xsite_table == 0x0 ) { 
	std::cerr << "Couldnt allocate site table\n" ;
	exit(1);
      }
    
      pad = 0;
      if( (ptrdiff_t)xsite_table % Cache::CacheLineSize != 0 ) { 
	pad = Cache::CacheLineSize - (ptrdiff_t)xsite_table % Cache::CacheLineSize;
      }
      site_table = (int *)((unsigned char *)xsite_table + pad);
      
      /* Loop through sites - you can choose your path below */
      /* The ordering below is checkerboarded in x, time fastest 
	 which should match exactly rb3 in QDP++ when compiled in cb3d mode */

      for(int p=0; p < 2; p++) { 
	for(int z=0; z < tot_size[2]; z++) {
	  for(int y=0; y < tot_size[1]; y++) { 
	    for(int x=0; x < tot_size[0]/2; x++) {
	      for(int t=0; t < tot_size[3]; t++) { 
		
		int coord[4];

		coord[0] = 2*x+p;
		coord[1] = y;
		coord[2] = z; 
		coord[3] = t;
	      
		/* Get the site and N-parity of the chosen victim */
		int  qdp_index = getLinearSiteIndex(coord); /* get the lexico index */
		int my_index = myLinearSiteIndex3D(coord);
		
		site_table[my_index]=qdp_index;

	      }
	    }
	  }
	}
      }
    
      /* Get the offsets needed for neighbour comm. */
      /* soffsets(position,direction,isign,cb)   */ 
      /*  where  isign    = +1 : plus direction */
      /*                  =  0 : negative direction */
      /*         cb       =  0 : even lattice (includes origin) */
      /*                  = +1 : odd lattice (does not include origin) */
      /* the offsets cotain the current site, i.e the neighbour for site i  */
      /* is  shift_table(i,dir,cb,mu) and NOT  i + soffset(..)    */
      
      /* Loop over directions and sites, building up shift tables */
      
      
      /* Loop through the sites linearly */
      for(int cb3=0; cb3 < 2; cb3++) { 
	for(int site=0; site < total_vol_cb; ++site) {
	  
	  int fcoord[4], bcoord[4];
	  int blinear, flinear;
	  int ipos;
	  
	  int my_index = cb3*total_vol_cb + site;
	  int qdp_index = site_table[ my_index ];
	  
	  /* Get the global site coords from the node and linear index */
	  int coord[4];
	  
	  getSiteCoords(coord, 0, qdp_index); 
	  for(dir=0; dir < Nd3; dir++) {	
	    
	    /* Backwards displacement*/
	    offs(bcoord, coord, dir, -1);
	    blinear = getLinearSiteIndex(bcoord);
	  
	    /* Forward displacement */
	    offs(fcoord, coord, dir, +1);
	    flinear = getLinearSiteIndex(fcoord);

	    
	    /* Gather */
	    shift_table[dir+Nd3*my_index ] = blinear; /* Linear index for psi, gauge */
	    shift_table[dir+Nd3*(my_index+total_vol)] = flinear; /* Linear index, psi or gauge */

	  }
	}
      }
    } 

  };

} // Namespace

#endif
