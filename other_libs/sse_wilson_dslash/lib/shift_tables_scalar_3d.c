/* $Id: shift_tables_scalar_3d.c,v 1.4 2008-03-05 19:45:12 bjoo Exp $

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

  static int Nd3 = 3;

  /* SIte table pointer */
  static int* xsite_table_3d;

  /* Aligned */
  int *site_table_3d;

  typedef struct { 
    int cb3;
    int linearcb3;
  } InvTab;


  /* Number of dimensions */
  static int getNumDim()
  {
    return (int)(4);
  }


/* Total problem size */
  static int tot_size_3d[4];
  static int total_vol_3d = -1;
  int total_vol_3d_cb = -1;

  static void setLattSize(const int size[])
  {
    int i;
    for(i=0; i < getNumDim(); ++i) {
      tot_size_3d[i] = size[i];
    }
    total_vol_3d = size[0];
    for(i=1; i < getNumDim(); ++i) {
      total_vol_3d *= size[i];
    }

    total_vol_3d_cb = total_vol_3d / 2;
  }

  static int* getLattSize()
  {
    return tot_size_3d;
  }

  static void crtesn3d(int ipos, const int latt_size[], int coord[] )
  {
  
    int Ndim=3;
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

  static int local_site(const coord[4], const latt_size[4])
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

  static int myLinearSiteIndex3D(const int gcoords[]) 
  {
    int mu;
    int subgrid_cb_nrow[4];
    int subgrid_cb_coord[4];
    int cb3;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = getLattSize()[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */

    cb3=0;
    for(mu=0; mu < Nd3; ++mu) { 
      cb3 += gcoords[mu];
    }
    cb3 &=1;
    
    subgrid_cb_coord[0] = (gcoords[0]/2)% subgrid_cb_nrow[0];
    for(mu=1; mu < 4; mu++) { 
      subgrid_cb_coord[mu] = gcoords[mu] % subgrid_cb_nrow[mu];
    }

    return local_site(subgrid_cb_coord, subgrid_cb_nrow) + cb3*total_vol_3d_cb;
  }

  // This is not needed as it can be done transitively:
  // ie lookup the QDP index and then lookup the coord with that 
  static void mySiteCoords3D(int gcoords[], int node, int linearsite)
  {
    int mu;
    int subgrid_cb_nrow[4];
    int tmp_coord[4];
    int cb3,cbb3;
    int my_node = QMP_get_node_number();

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = getLattSize()[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */

    /* Single processor -- all coords 0 */
    for(mu=0; mu < 4; mu++) { 
      gcoords[mu] = 0;
    }
    
    cb3=linearsite/total_vol_3d_cb;

    crtesn3d(linearsite % total_vol_3d_cb, subgrid_cb_nrow, tmp_coord);

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

  /* Offset by 1 in direction dir */
  static void offs(int temp[], const int coord[], int mu, int isign)
  {
    int i;
    
    for(i=0; i < getNumDim(); ++i)
      temp[i] = coord[i];
    
    /* translate address to neighbour */
    temp[mu] = (temp[mu] + isign + 2*getLattSize()[mu]) % getLattSize()[mu];
  }
  
  /* This is an Nd parity */
  static int parity(const int coord[])
  {
    int m;
    int sum = 0;
    
    for(m=0; m < Nd3; ++m)
      sum += coord[m];
    
    return sum % 2;
  }


 
  int* make_shift_tables_3d(const int nrow[],
			    void (*QDP_getSiteCoords)(int coord[], int node, int linearsite),
			    
			    int (*QDP_getLinearSiteIndex)(const int coord[]))
  { 
    int dir; 
    int coord[4];
    int linear;
    int backward=0;
    int forward =1;
    int *shift_table;
    
    int x,y,z,t;
    int cb3;
    
    int site;
    int p;
    int qdp_index; 
    int my_index;
    
    /* Set the lattice size, get total volume and checkerboarded volume */
    setLattSize(nrow);
    

    /* Allocate the shift table */
    /* Nd3 directions
       total vol sites
       2 types (FWD, BKWD)
    */
    if ((shift_table = (int *)malloc(Nd3*total_vol_3d*2*sizeof(int))) == 0) {
      fprintf(stderr,"init_sse_su3dslash: could not initialize shift_table\n");
      exit(1);
    }
    
    /* Allocate the site table and the shift table */
    /* Now I want to build the site table */
    /* I want it cache line aligned? */
    xsite_table_3d = (int *)malloc(sizeof(int)*total_vol_3d+63);
    if(xsite_table_3d == 0x0 ) { 
      fprintf(stderr,"Couldnt allocate site table\n");
      exit(1);
    }
    site_table_3d = (int *)((((ptrdiff_t)(xsite_table_3d))+63L)&(-64L));
    
    
    /* Loop through sites - you can choose your path below */
    /* The ordering below is checkerboarded in x, time fastest 
       which should match exactly rb3 in QDP++ when compiled in cb3d mode */

    for(p=0; p < 2; p++) { 
      for(z=0; z < nrow[2]; z++) {
	for(y=0; y < nrow[1]; y++) { 
	  for(x=0; x < nrow[0]/2; x++) {
	    for(t=0; t < nrow[3]; t++) { 
	    
	      coord[0] = 2*x+p;
	      coord[1] = y;
	      coord[2] = z; 
	      coord[3] = t;
	      
	      /* Get the site and N-parity of the chosen victim */
	      qdp_index = QDP_getLinearSiteIndex(coord); /* get the lexico index */
	      my_index = myLinearSiteIndex3D(coord);

	      site_table_3d[my_index]=qdp_index;

	     
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
    for(cb3=0; cb3 < 2; cb3++) { 
      for(site=0; site < total_vol_3d_cb; ++site) {
	
	int fcoord[4], bcoord[4];
	int blinear, flinear;
	int ipos;
	
	my_index = cb3*total_vol_3d_cb + site;
	qdp_index = site_table_3d[ my_index ];
	
	/* Get the global site coords from the node and linear index */
	QDP_getSiteCoords(coord, 0, qdp_index); 
	for(dir=0; dir < Nd3; dir++) {	
	  
	  /* Backwards displacement*/
	  offs(bcoord, coord, dir, -1);
	  
	  blinear = QDP_getLinearSiteIndex(bcoord);
	  
	  /* Forward displacement */
	  offs(fcoord, coord, dir, +1);
	  flinear = QDP_getLinearSiteIndex(fcoord);
	  
	  /* Gather */
	  shift_table[dir+Nd3*my_index ] = blinear; /* Linear index for psi, gauge */
	  shift_table[dir+Nd3*(my_index+total_vol_3d)] = flinear; /* Linear index, psi or gauge */
	}
      }
    }
    
    
    return shift_table;
  } 

int forward_neighbor_3d(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + Nd3*(mysite + total_vol_3d)];
}

int backward_neighbor_3d(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + Nd3*mysite];
}

void free_shift_tables_3d(int **table) {
  free(*table);
  free(xsite_table_3d);
}


#ifdef __cplusplus
}
#endif
