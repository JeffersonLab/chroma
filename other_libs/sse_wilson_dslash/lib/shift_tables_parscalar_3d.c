

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "qmp.h"
#include "shift_tables_parscalar.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Unaligned Offset Table */
  static int* xoffset_table_body_3d;

  /* EXPORTED: Aligned Offset Table */
  int* offset_table_body_3d;

  /* A struct for the inverse site table. This may be eliminated soon */
  typedef struct { 
    int cb3;
    int linearcb3;
  } InvTab;

  /* Unaligned Site Table */
  static int* xsite_table_3d;

  /* EXPORTED: Aligned Site Table */
  /* This runs as site_table_3d[index] = qdp_site_index */
  /* Now 'index' runs in checkerboarded order           */
  /* qdp_site_index is the offset into QDP arrays for index */
     
  int* site_table_3d;


  /* Total Subgrid Volume on Core */
  /* EXPORTED:                    */
  int subgrid_vol_3d = -1;

  /* CB Volume on Core */
  /* EXPORTED */
  int subgrid_vol_cb_3d = -1;

  /* The number of dimensions once and for all */
  /* EXPORTED */
  int Nd3 = 3;

  /* Number of dimensions */
  static int getNumDim()
  {
    return (int)(QMP_get_logical_number_of_dimensions());
  }


  /* Subgrid lattice size */
  static int* getSubgridSize()
  {
    static int first = 1;
    static int subgrid_size[4];

    if (first == 1) {
      int i;
      for(i=0; i < getNumDim(); ++i) {
	subgrid_size[i] = QMP_get_subgrid_dimensions()[i];
      }

      /* Why do we multiply by 2 after QMP_get_subgrid_dimensions() */
      subgrid_size[0]  *= 2;
      
      first = 0;
    }
  
    return subgrid_size;
  }


  /* Some functions */
  /* 
     int myLinearSiteIndex(int gcoords[], int *node)
     
     Maps globalCoordinates into a local site index.
     This site index can be used to directly index external half spinors
     and is the equivalent of getLinearSiteIndex for my internal layout.
     Valid only on processor where the coords live.
  */

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


  static int local_site(const int coord[4], const int latt_size[4])
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
      subgrid_cb_nrow[mu] = getSubgridSize()[mu];
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

    return local_site(subgrid_cb_coord, subgrid_cb_nrow) + cb3*subgrid_vol_cb_3d;
  }


  // This is not needed as it can be done transitively:
  // ie lookup the QDP index and then lookup the coord with that 
  static void mySiteCoords3D(int gcoords[], int node, int linearsite)
  {
    int mu;
    int subgrid_cb_nrow[4];
    int tmp_coord[4];
    int cb3,cbb3;
    int* log_coords=QMP_get_logical_coordinates_from(node);
    int my_node = QMP_get_node_number();

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = getSubgridSize()[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */


    for(mu=0; mu < 4; mu++) { 
      gcoords[mu] = log_coords[mu]*getSubgridSize()[mu];

    }
    
    cb3=linearsite/subgrid_vol_cb_3d;

    crtesn3d(linearsite % subgrid_vol_cb_3d, subgrid_cb_nrow, tmp_coord);

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





  /* Total problem size */
  /* File scope so no conflict with 4D */
  static int* getLattSize()
  {
    static int first = 1;
    static int tot_size[4];

    if (first == 1) {
      
      const int* phys_size = QMP_get_logical_dimensions();
      int i;
      
      for(i=0; i < getNumDim(); ++i) {
	tot_size[i] = getSubgridSize()[i]*phys_size[i];
      }
      
      first = 0;
    }
    
    return tot_size;
  }



  /* Offset by 1 in direction dir */
  static void offs(int temp[], const int coord[], int mu, int isign)
  {
    int i;
    
    for(i=0; i < getNumDim(); ++i) {
      temp[i] = coord[i];
    }
    /* translate address to neighbour */
    temp[mu] = (temp[mu] + isign + 2*getLattSize()[mu]) % getLattSize()[mu];
  }


  /* This is now the Nd3 parity */
  static int parity(const int coord[])
  {
    int m;
    int sum = 0;
    
    for(m=0; m < Nd3; ++m)
      sum += coord[m];
    
    return sum & 1;
  }


 

  /* This one does al the hard work */
  void make_shift_tables_3d(int bound[2][4][3],
			    void (*QDP_getSiteCoords)(int coord[], int node, int linearsite), 
			    int (*QDP_getLinearSiteIndex)(const int coord[]),

			    int (*QDP_getNodeNumber)(const int coord[]))
  { 
    volatile int dir,i; 
    const int my_node = QMP_get_node_number();
    
    int coord[4];
    int gcoord[4];
    int gcoord2[4];

    int linear;
    int **shift_table;
    int x,y,z,t;
    int *subgrid_size = getSubgridSize();
    int mu;
    
    int cb3;
    const int *node_coord  = QMP_get_logical_coordinates();
    int p;
    int site,index;

    InvTab *xinvtab;
    InvTab *invtab;

    int qdp_index;
    int my_index;

    /* Setup the subgrid volume for ever after */
    /* Nb getNumDim() here returns the number of QMP dimensions */
    /* This is potentially undesirable - maybe I should just hardwire this to be 4 */
    /* This sets up a global, allowing subgrid_vol_cb3 to be used ever after -- Yuck! */
    subgrid_vol_3d = 1;
    for(i=0; i < getNumDim(); ++i) {
      subgrid_vol_3d *= getSubgridSize()[i]; 
    }
    

    /* Get the checkerboard size for ever after */
    /* Again this sets up a global. Yuck */
    subgrid_vol_cb_3d = subgrid_vol_3d / 2;


    /* Now I want to build the site table */
    /* I want it cache line aligned? */
    xsite_table_3d = (int *)malloc(sizeof(int)*subgrid_vol_3d+63);
    if(xsite_table_3d == 0x0 ) { 
      QMP_error("Couldnt allocate site table");
      QMP_abort(1);
    }
    site_table_3d = (int *)((((ptrdiff_t)(xsite_table_3d))+63L)&(-64L));

    /* This is an inverse table */
    xinvtab = (InvTab *)malloc(sizeof(InvTab)*subgrid_vol_3d+63);
    if(xinvtab == 0x0 ) { 
      QMP_error("Couldnt allocate site table");
      QMP_abort(1);
    }
    invtab = (InvTab *)((((ptrdiff_t)(xinvtab))+63L)&(-64L));

    /* Inversity of functions check:
       Check that myLinearSiteIndex3D is in fact the inverse
       of mySiteCoords3D, and that QDP_getSiteCoords is the
       inverse of QDP_linearSiteIndex()
    */
    for(p=0; p < 2; p++) {
      for(site=0; site < subgrid_vol_cb_3d; site++) { 
	
	/* Linear site index */
	my_index = site + subgrid_vol_cb_3d*p;
	QDP_getSiteCoords(gcoord, my_node, my_index);
	linear=QDP_getLinearSiteIndex(gcoord);

	if( linear != my_index ) { 
	  printf("P%d cb=%d site=%d : QDP_getSiteCoords not inverse of QDP_getLinearSiteIndex(): my_index=%d linear=%d\n", my_node, p,site, my_index,linear);
	}

	mySiteCoords3D(gcoord, my_node, my_index);
	linear=myLinearSiteIndex3D(gcoord);

	if( linear != my_index ) { 
	  printf("P%d cb=%d site=%d : mySiteCoords3D not inverse of myLinearSiteIndex3D(): my_index=%d linear=%d\n", my_node, p,site, my_index,linear);
	}
      }
    }


    /* Loop through the sites in some order.
       This is an ordering that should coincide with QDP++'s rb3
       subset, when QDP++ is compiled in cb3d order */

    for(p=0; p < 2; p++) { 	      
      for(z=0; z < subgrid_size[2]; z++) { 
	for(y=0; y < subgrid_size[1]; y++) { 
	  for(x=0; x < subgrid_size[0]/2; x++) { 
	    for(t=0; t < subgrid_size[3]; t++) { 
	      
	      coord[0] = 2*x + p;	      
	      coord[1] = y;
	      coord[2] = z;
	      coord[3] = t;

	      /* Make global */
	      for(i=0; i < 4; i++) { 
		coord[i] += subgrid_size[i]*node_coord[i];
	      }

	      /* Both these indices serve as an index into something 
		 of lattice size. */

	      /* Index of coordinate -- NB this is not lexicographic
		 but takes into account checkerboarding in QDP++ */
	      qdp_index = QDP_getLinearSiteIndex(coord);

	      /* Index of coordinate in my layout. -- NB this is not lexicographic
		 but takes into account my 3D checkerbaording */
	      my_index = myLinearSiteIndex3D(coord);
	      

	      site_table_3d[my_index] = qdp_index;

	      cb3=parity(coord);
	      linear = my_index%subgrid_vol_cb_3d;

	      invtab[qdp_index].cb3=cb3;
	      invtab[qdp_index].linearcb3=linear;
	    }
	  }
	}
      }
    }

    /* Site table transitivity check: 
       for each site, convert to index in cb3d, convert to qdp index
                convert qdp_index to coordinate
		convert coordinate to back index in cb3d
       Check that your cb3d at the end is the same as you 
       started with */
    for(p=0; p < 2; p++) { 
      for(site=0; site < subgrid_vol_cb_3d; site++) {

	/* My local index */
	my_index = site + subgrid_vol_cb_3d*p;
	
	/* Convert to QDP index */
	qdp_index = site_table_3d[ my_index ];

	/* Switch QDP index to coordinates */
	QDP_getSiteCoords(gcoord, my_node,qdp_index);

	/* Convert back to cb3d index */
	linear = myLinearSiteIndex3D(gcoord);

	/* Check new cb3d index matches the old cb3d index */
	if (linear != my_index) { 
	  printf("P%d The Circle is broken. My index=%d qdp_index=%d coords=%d,%d,%d,%d linear(=my_index?)=%d\n", my_node, my_index, qdp_index, gcoord[0],gcoord[1],gcoord[2],gcoord[3],linear);
	}
      }
    }


    /* Consistency check 2: Test mySiteCoords 3D
       for all 3d cb,cb3index convert to
                   cb3d linear index (my_index)
        convert to qdp_index (lookup in site table)

       Now convert qdp_index and my_index both to 
       coordinates. They should produce the same coordinates
    */
    for(p=0; p < 2; p++) { 
      for(site=0; site < subgrid_vol_cb_3d; site++) {

	/* My local index */
	my_index = site + subgrid_vol_cb_3d*p;
	mySiteCoords3D(gcoord, my_node, my_index);

	qdp_index = site_table_3d[ my_index ];
	QDP_getSiteCoords(gcoord2, my_node,qdp_index);
	
	for(mu=0 ; mu < 4; mu++) { 
	  if( gcoord2[mu] != gcoord[mu] ) {
	    printf("P%d: my_index=%d qdp_index=%d mySiteCoords=(%d,%d,%d,%d) QDPsiteCoords=(%d,%d,%d,%d)\n", my_node, my_index, qdp_index, gcoord[0], gcoord[1], gcoord[2], gcoord[3], gcoord2[0], gcoord2[1], gcoord2[2], gcoord2[3]);
	    continue;
	  }
	}
	
      }
    }

    /* Allocate the shift table */
    /* The structure is as follows: There are 4 shift tables in order:
       
       [ Table 1 | Table 2 | Table 3 | Table 4 ]
       Table 1: decomp_scatter_index[mu][site]
       Table 2: decomp_hvv_scatter_index[mu][site]
       Table 3: recons_mvv_gather_index[mu][site]
       Table 4: recons_gather_index[mu][site]
       
       Each table is indexed as:
       table[ type ][ dir-site ]  
       
       where dir-site is a combination of
       direction and site indices with direction running fastest.
       
       dir_site is computed as dir + Nd3*site 
       
       site runs as linear+subgrid_vol_cb*cb3 - although 
       this is also rolled into a linear loop over all the sites.
       
       The values stored in the index table:
       i) If the appropriate neighbour is on-site, we give the 
       linear part of its (cb3,linear) coordinate. The CB3 is 
       implied (essentially its a neigbour so it is opposite ours)
       or it is a received buffer, and may be the same as ours.
       In any case, the stored index will only index something
       on one checkerboard and so we don't need to supply the 
       checkerboard info.
       
       ii) If the neighbour is off site, we index into a tail.
       
       Our temporaries will be arranged as:
       
       [  body half spinors ][ Tail 1 half spinors ][ Tail 2 half spinors ] 
       Tail 1 gets sent backward and tail 2 receives from forward 
       or vice versa depending on parity. In any case DECOMPs always go to
       Tail 1, and RECONS-s always go to Tail 2.
       
       We count the sizes of the tails as we go along, and fill out the 
       bound array. Typically tail 1 indices look like subgrid_vol_cb+x
       While tail 2 indices look like 2*subgrid_vol_cb + x
    */
    
    /* 4 for the four types above: */
    if ((shift_table = (int **)malloc(4*sizeof(int*))) == 0 ) {
      QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
      QMP_abort(1);
      
    }
    
    /* Now the table for each of the 4 types */
    for(i=0; i < 4; i++) { 
      /* Nd3 for the  directions */
      if ((shift_table[i] = (int *)malloc(Nd3*subgrid_vol_3d*sizeof(int))) == 0) {
	QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
	QMP_abort(1);
      }
    }
    

    /* Initialize the boundary counters */
    for(cb3=0; cb3 < 2; cb3++) {
      for(dir=0; dir < Nd3 ; dir++) {
	bound[cb3][0][dir] = 0;	
	bound[cb3][1][dir] = 0;	
	bound[cb3][2][dir] = 0;	
	bound[cb3][3][dir] = 0;	
      }
    }



    /* Loop over All sites */
    for(cb3=0; cb3 < 2; cb3++) 
      for(site=0; site < subgrid_vol_cb_3d; ++site) {
	
	index = cb3*subgrid_vol_cb_3d + site;
	/* Get the QDP offset corresponding to my (cb,site) */
	qdp_index = site_table_3d[index]; 
	
	
	/* Get the global site coords from the node and linear index */
	QDP_getSiteCoords(coord, my_node, qdp_index);
       
	
	/* Loop through the directions */
	for(dir=0; dir < Nd3; dir++) {
	  
	  int fcoord[4], bcoord[4];
	  int fnode, bnode;
	  int blinear, flinear;
	  
	  /* Find forward and backward neighbours. Get both 
	     the node info and the linear coordinates. 
	     NB The linear coordinates are the absolute 
	     lexico ones. */
	  
	  /* Backwards displaced coordinate & Node */
	  offs(bcoord, coord, dir, -1);
	  bnode   = QDP_getNodeNumber(bcoord);  /* Its node */
	  blinear = QDP_getLinearSiteIndex(bcoord); /* Its QDP++ linear offset */
	  
	  /* Forward displaced coordinate & Node */
	  offs(fcoord, coord, dir, +1);
	  fnode   = QDP_getNodeNumber(fcoord);
	  flinear = QDP_getLinearSiteIndex(fcoord);
	  
	  
	  /* Scatter:  decomp_{plus,minus} */
	  /* Operation: a^F(shift(x,type=0),dir) <- decomp(psi(x),dir) */ 
	  /* Send backwards - also called a receive from forward */
	  
	if (bnode != my_node) {      
	  /* Offnode */
	  /* Append to Tail 1, increase boundary count */
	  shift_table[DECOMP_SCATTER][dir+Nd3*index] 
	    = subgrid_vol_cb_3d + bound[1-cb3][DECOMP_SCATTER][dir];
	  
	  bound[1-cb3][DECOMP_SCATTER][dir]++;
	  
	}
	else {                                           
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup */
	  shift_table[DECOMP_SCATTER][dir+Nd3*index] =
	    invtab[blinear].linearcb3; // My linear (cb implied)
	}
	
	/* Scatter:  decomp_hvv_{plus,minus} */
	/* Operation:  a^B(shift(x,type=1),dir) <- U^dag(x,dir)*decomp(psi(x),dir) */
	/* Send forwards - also called a receive from backward */
	if (fnode != my_node) {
	  /* Offnode */
	  /* Append to Tail 1, increase boundary count */
	  shift_table[DECOMP_HVV_SCATTER][dir+Nd3*index]           
	    = subgrid_vol_cb_3d + bound[1-cb3][DECOMP_HVV_SCATTER][dir];
	  
	  bound[1-cb3][DECOMP_HVV_SCATTER][dir]++;                  
	  
	}
	else {
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup */
	  shift_table[DECOMP_HVV_SCATTER][dir+Nd3*index] =           /* Onnode */
	    invtab[flinear].linearcb3;
	}
	
	/* Gather:  mvv_recons_{plus,minus} */
	/* Operation:  chi(x) <-  \sum_dir U(x,dir)*a^F(shift(x,type=2),dir) */
	/* Receive from forward */
	if (fnode != my_node) {
	  /* Offnode */
	  /* Append to Tail 2, increase boundary count */
	  
	  shift_table[RECONS_MVV_GATHER][dir+Nd3*index] =
	    2*subgrid_vol_cb_3d + (bound[cb3][RECONS_MVV_GATHER][dir]);
	  
	  bound[cb3][RECONS_MVV_GATHER][dir]++;
	  
	}
	else {
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup. Note this is a recons post shift,
	     so the linear coordinate to invert is mine rather than the neighbours */
	  shift_table[RECONS_MVV_GATHER][dir+Nd3*index] =
	    invtab[qdp_index].linearcb3;
	}
	
	
	/* Gather:  recons_{plus,minus} */
	/* Operation:  chi(x) +=  \sum_dir recons(a^B(shift(x,type=3),dir),dir) */
	/* Receive from backward */
	if (bnode != my_node) {
	  
	  shift_table[RECONS_GATHER][dir+Nd3*index] = 
	    2*subgrid_vol_cb_3d + bound[cb3][RECONS_GATHER][dir];
	  
	  bound[cb3][RECONS_GATHER][dir]++;
	}
	else {
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup. Note this is a recons post shift,
	     so the linear coordinate to invert is mine rather than the neighbours */
	  
	  shift_table[RECONS_GATHER][dir+Nd3*index] =
	    invtab[qdp_index].linearcb3;
	}
	
      } 
    }
  

    
    /* Now I want to make the offset table into the half spinor temporaries */
    /* The half spinor temporaries will look like this:
       
       dir=0 [ Body Half Spinors ][ Tail 1 Half Spinors ][ Tail 2 Half Spinors ]
       dir=1 [ Body Half Spinors ][ Tail 1 Half Spinors ][ Tail 2 Half Spinors ]
       ...
       
       And each of these blocks of half spinors will be sized to vol_cb
       sites (ie half volume only).  The shift_table() for a given site and
       direction indexes into one of these lines. So the offset table essentially
       delineates which line one picks, by adding an offset of 
       3*subgrid_vol_cb*dir 
       To the shift. The result from offset table, can be used directly as a
       pointer displacement on the temporaries.
       
    */
    /* The 4 is for the 4 types, Nd 3 is for the 3 dirs, subgrid_vol_3d 
       is for all the sites */
    xoffset_table_body_3d = (int *)malloc(Nd3*4*subgrid_vol_3d*sizeof(int)+63);
    if( xoffset_table_body_3d == 0 ) {
      QMP_error("init_wnxtsu3dslash: could not initialize offset_table[i]");
      QMP_abort(1);
    }
    /* This is the bit what aligns straight from AMD Manual */
    offset_table_body_3d = (int *)((((ptrdiff_t)(xoffset_table_body_3d)) + 63L) & (-64L));
    
    /* Loop over all the sites */
    for(site=0; site < subgrid_vol_3d; site++) { 
      for(dir=0; dir < Nd3; dir++) { 
	
	/* Here site is a local site */
	offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*DECOMP_SCATTER ) ] = 
	  shift_table[DECOMP_SCATTER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
	
	offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*DECOMP_HVV_SCATTER ) ] = 
	  shift_table[DECOMP_HVV_SCATTER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
	
	offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*RECONS_MVV_GATHER ) ] = 
	  shift_table[RECONS_MVV_GATHER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
	
	offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*RECONS_GATHER ) ] = 
	  shift_table[RECONS_GATHER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
	
      }
    }
    
    /* Free shift table - it is no longer needed. We deal solely with offsets */
    for(i=0; i < 4; i++) { 
      free( (shift_table[i]) );
    }
    free( shift_table );
    
    /* Free the inverse site table, it is no longer needed */
    free( xinvtab );
    
  } 
  

void free_shift_tables_3d(void) 
{
  /* Free the offset and site tables - free their aligned variants */
  free( xoffset_table_body_3d );
  free( xsite_table_3d );

}

#ifdef __cplusplus
}
#endif
