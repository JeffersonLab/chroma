/* $Id: shift_tables_scalar.c,v 1.5 2008-03-05 19:45:12 bjoo Exp $

/* Set the offset tables used by the 32-bit and 64-bit single node dslash */

/* make_shift_tables(latt_size) */
/* Get the offsets needed for neighbour comm. */
/* soffsets(position,direction,isign,cb)   */ 
/*  where  isign    = +1 : plus direction */
/*                  =  0 : negative direction */
/*         cb       =  0 : even lattice (includes origin) */
/*                  = +1 : odd lattice (does not include origin) */
/* the offsets cotain the current site, i.e the neighbour for site i  */
/* is  soffsets(i,dir,cb,mu) and NOT  i + soffset(..)    */
  
#include <shift_tables_scalar.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


  static int Nd=4;
  static int* xsite_table;
  int* site_table;


  /* Total problem size */
  static int tot_size[4];
  static int total_vol = -1;
  int total_vol_cb = -1;

  static void setLattSize(const int size[])
  {
    int i;
    for(i=0; i < 4; ++i) {
      tot_size[i] = size[i];
    }
    total_vol = size[0];
    for(i=1; i < 4; ++i) {
      total_vol *= size[i];
    }

    total_vol_cb = total_vol / 2;
  }


  static int* getLattSize()
  {
    return tot_size;
  }


  /* Decompose lexicographic site ipos into lattice coordinates */
  static void crtesn4d(int ipos, const int latt_size[], int coord[] )
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

  /* Calculates the lexicographic site index from the coordinate of a site */
  static int local_site4d(int coord[], int latt_size[])
  {
    int order = 0;
    int mmu;
    
    for(mmu=4-1; mmu >= 1; --mmu) {
      order = latt_size[mmu-1]*(coord[mmu] + order);
    }
    order += coord[0];
    
    return order;
  }

  static int myLinearSiteIndex4D(const int gcoords[]) 
  {
    int mu;
    int subgrid_cb_nrow[4];
    int subgrid_cb_coord[4];
    int cb;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = getLattSize()[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */

    cb=0;
    for(mu=0; mu < Nd; ++mu) { 
      cb += gcoords[mu];
    }
    cb &=1;
    
    subgrid_cb_coord[0] = (gcoords[0]/2)% subgrid_cb_nrow[0];
    for(mu=1; mu < 4; mu++) { 
      subgrid_cb_coord[mu] = gcoords[mu] % subgrid_cb_nrow[mu];
    }

    return local_site4d(subgrid_cb_coord, subgrid_cb_nrow) + cb*total_vol_cb;
  }

  // This is not needed as it can be done transitively:
  // ie lookup the QDP index and then lookup the coord with that 
  static void mySiteCoords4D(int gcoords[], int node, int linearsite)
  {
    int mu;
    int subgrid_cb_nrow[4];
    int tmp_coord[4];
    int cb,cbb;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = getLattSize()[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */

    /* Base coordinate single processor: 0,0,0,0 always */
    for(mu=0; mu < 4; mu++) { 
      gcoords[mu] = 0;
    }
    
    cb=linearsite/total_vol_cb;

    crtesn4d(linearsite % total_vol_cb, subgrid_cb_nrow, tmp_coord);

    // Add on position within the node
    // NOTE: the cb for the x-coord is not yet determined
    gcoords[0] += 2*tmp_coord[0];
    for(mu=1; mu < 4; ++mu) {
      gcoords[mu] += tmp_coord[mu];
    }

    cbb = cb;
    for(mu=1; mu < 4; ++mu) {
      cbb += gcoords[mu];
    }
    gcoords[0] += (cbb & 1);
  }

 

  /* Offset by 1 in direction dir */
  static void offs(int temp[], const int coord[], int mu, int isign)
  {
    int i;
    
    for(i=0; i < 4; ++i) {
      temp[i] = coord[i];
    }

    /* translate address to neighbour */
    temp[mu] = (temp[mu] + isign + 2*getLattSize()[mu]) % getLattSize()[mu];
  }

  /* This is a 4D parity */
  static int parity(const int coord[])
  {
    int m;
    int sum = 0;
    
    for(m=0; m < 4; ++m) {
      sum += coord[m];
    }
    return sum % 2;
  }


 
  int* make_shift_tables(const int nrow[],
			 void (*QDP_getSiteCoords)(int coord[], int node, int linearsite),
			 
			 int (*QDP_getLinearSiteIndex)(const int coord[]))
  
{ 
  int dir; 
  int coord[4];
  int linear;
  int Nd = 4;
  int backward=0;
  int forward =1;
  int *shift_table;

  int x,y,z,t;
  int cb;
  int site;

  int qdp_index; 
  int my_index;
  int p;

  /* Set the lattice size, get total volume and checkerboarded volume */
  setLattSize(nrow);

  /* Determine what is the starting site for each color (checkerboard) */
  // icolor_start[0] = 0;
  //icolor_start[1] = total_vol_cb;

  /* Allocate the shift table */
  if ((shift_table = (int *)malloc(4*total_vol*2*sizeof(int))) == 0) {
    fprintf(stderr,"init_sse_su3dslash: could not initialize shift_table\n");
    exit(1);
  }

  /* Allocate the site table and the shift table */
  /* Now I want to build the site table */
  /* I want it cache line aligned? */
  xsite_table = (int *)malloc(sizeof(int)*total_vol+63);
  if(xsite_table == 0x0 ) { 
    fprintf(stderr,"Couldnt allocate site table\n");
    exit(1);
  }
  site_table = (int *)((((ptrdiff_t)(xsite_table))+63L)&(-64L));

  /* Loop through sites - you can choose your path below */
  /* This is the ordering for QDP++ in CB2 mode. - In this mode
     the site tables of QDP rb2 subset and mine should match */
  
  for(p=0; p < 2; p++) { 	    
    for(t=0; t < nrow[3]; t++) { 
      for(z=0; z < nrow[2]; z++) {
	for(y=0; y < nrow[1]; y++) {     
	  for(x=0; x < nrow[0]/2; x++) {

	    coord[0] = 2*x+p;
	    coord[1] = y;
	    coord[2] = z; 
	    coord[3] = t;
	    
	    /* Get the site and N-parity of the chosen victim */
	    qdp_index = QDP_getLinearSiteIndex(coord); /* get the lexico index */
	    my_index = myLinearSiteIndex4D(coord);
	    
	    /* Add lexico site into site_table, for current cb3 and linear */
	    /* Map (cb3, linear) -> lexico */
	    site_table[ my_index ] = qdp_index;

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
  for(cb=0; cb < 2; cb++) {
    for(site = 0; site < total_vol_cb; ++site) { 
      int fcoord[4], bcoord[4];
      int blinear, flinear;
      
      my_index = cb*total_vol_cb + site;
      
      qdp_index = site_table[ my_index ];

      QDP_getSiteCoords(coord, 0, qdp_index); 

      for(dir=0; dir < 4; dir++) {

	/* Backwards displacement*/
	offs(bcoord, coord, dir, -1);
	blinear = QDP_getLinearSiteIndex(bcoord);

	/* Forward displacement */
	offs(fcoord, coord, dir, +1);
	flinear = QDP_getLinearSiteIndex(fcoord);


	/* Gather */
	shift_table[dir+Nd*my_index ] = blinear;
	shift_table[dir+Nd*(my_index+total_vol)] = flinear;
      }
    }
  }

  return shift_table;
} 

int forward_neighbor(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + 4*(mysite + total_vol)];
}

int backward_neighbor(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + 4*mysite];
}

void free_shift_tables(int **table) {
  free(*table);
  free(xsite_table);
}

/* Total volume */
int getSubgridVol()
{
  return total_vol;
}

/* Total CB volume */
int getSubgridVolCB()
{
  return total_vol_cb;
}

#ifdef __cplusplus
}
#endif
