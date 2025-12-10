/* $Id: shift_tables_parscalar.c,v 1.12 2008-07-30 20:25:04 bjoo Exp $ */


/* both of these must be called before the P4 dslash is called */
/* computes the scatter/gather table for the P4 32-bit parallel dslash */

/* make_shift_tables(
 * int *shift,               the external table to write the new shift table to
 * int subgrid_vol_cb         # of sites that a call to the dslash will run over per node
 * int nrow[]                original problem size (before checkerboarding)
 * int subgrid_cb_nrow[]     number of sites per subgrid per cb for each direction
 * int bound[]    array of size Nd that holds the number of boundaries per direction
 * int Nd         number of directions
 * )
 * the new shift table will be accessible with hte following index order:
 * shift(
      dir   --> direction
      site  --> site
      forward/backward --> 0 backward, 1 forward
      cb   --> cb
      type --> 0: backward=scatter, forward=gather  1: backward=gather, forward=scatter)

 * decomp uses back scatter type = 0 forw/back = 0
 * decomp_hvv uses forward scatter type = 1 forw/back = 1
 * mvv_recons uses forward gather type = 0 forw/back = 1
 * recons uses back gather type = 1 forw/back = 0

 next version needs to reorder these indices in this file and then propagate to make cb slowest varying and to ditch this forward backward and type stuff (which makes mathematical sense) and go for shift table 0, 1, 2, 3....that might help a little ..... 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "qmp.h"
#include "shift_tables_parscalar.h"

#ifdef __cplusplus
extern "C" {
#endif



  /* Offset table into the half spinor array */
  /* Unaligned */
  static halfspinor_array** xoffset_table;

  /* Aligned version */
  halfspinor_array** offset_table;

  /* A struct for the inverse of the shift table - may be eliminated in 
     future */
  typedef struct { 
    int cb;
    int linearcb;
  } InvTab4;

  /* Unaligned Site Table */
  static int* xsite_table;
  
  /* EXPORTED: Aligned Site Table */
  int* site_table;


  /* Max machine size */
  int subgrid_vol = -1;
  int subgrid_vol_cb = -1;

  int Nd = 4;


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

      
      subgrid_size[0] *= 2;
      
      first = 0;
    }
  
    return subgrid_size;
  }


  /* Total problem size */
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

  static int local_site4d(const coord[4], const latt_size[4])
  {
    int order = 0;
    int mmu;
    
    // In the 4D Case: t+Lt(x + Lx(y + Ly*z)
    // essentially  starting from i = dim[Nd-2]
    //  order =  latt_size[i-1]*(coord[i])
    //   and need to wrap i-1 around to Nd-1 when it gets below 0
    for(mmu=3; mmu >= 1; --mmu) {
      order = latt_size[mmu-1]*(coord[mmu] + order);
    }
    
    order += coord[ 0 ]; /* X is fastest running */
    
    return order;
  }

  static int myLinearSiteIndex4D(const int gcoords[]) 
  {
    int mu;
    int subgrid_cb_nrow[4];
    int subgrid_cb_coord[4];
    int cb;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = getSubgridSize()[mu];
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

    return local_site4d(subgrid_cb_coord, subgrid_cb_nrow) + cb*subgrid_vol_cb;
  }

  // This is not needed as it can be done transitively:
  // ie lookup the QDP index and then lookup the coord with that 
  static void mySiteCoords4D(int gcoords[], int node, int linearsite)
  {
    int mu;
    int subgrid_cb_nrow[4];
    int tmp_coord[4];
    int cb,cbb;
    int* log_coords=QMP_get_logical_coordinates_from(node);
    int my_node = QMP_get_node_number();

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = getSubgridSize()[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */


    for(mu=0; mu < 4; mu++) { 
      gcoords[mu] = log_coords[mu]*getSubgridSize()[mu];

    }
    
    cb=linearsite/subgrid_vol_cb;

    crtesn4d(linearsite % subgrid_vol_cb, subgrid_cb_nrow, tmp_coord);

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
    
    for(i=0; i < getNumDim(); ++i)
      temp[i] = coord[i];
    
    /* translate address to neighbour */
    temp[mu] = (temp[mu] + isign + 2*getLattSize()[mu]) % getLattSize()[mu];
  }
  
  
  static int parity(const int coord[])
  {
    int m;
    int sum = 0;
    
    for(m=0; m < 4; ++m)
      sum += coord[m];
    
    return sum & 1;
  }
  
  
  void make_shift_tables(int bound[2][4][4], halfspinor_array* chi1,
			 halfspinor_array* chi2,

			 halfspinor_array* recv_bufs[2][4],
			 halfspinor_array* send_bufs[2][4],

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
  int offset;

  int cb;
  const int *node_coord  = QMP_get_logical_coordinates();
  int p;
  int site, index;

  InvTab4 *xinvtab;
  InvTab4 *invtab;

  int qdp_index;
  int my_index;
  int num;
  int offsite_found;

  /* Setup the subgrid volume for ever after */
  subgrid_vol = 1;
  for(i=0; i < getNumDim(); ++i) {
    subgrid_vol *= getSubgridSize()[i]; 
  }
  
  /* Get the checkerboard size for ever after */
  subgrid_vol_cb = subgrid_vol / 2;

  /* Now I want to build the site table */
  /* I want it cache line aligned? */
  xsite_table = (int *)malloc(sizeof(int)*subgrid_vol+63L);
  if(xsite_table == 0x0 ) { 
    QMP_error("Couldnt allocate site table");
    QMP_abort(1);
  }

  site_table = (int *)((((ptrdiff_t)(xsite_table))+63L)&(-64L));

  xinvtab = (InvTab4 *)malloc(sizeof(InvTab4)*subgrid_vol + 63L);
  if(xinvtab == 0x0 ) { 
    QMP_error("Couldnt allocate site table");
    QMP_abort(1);
  }
  invtab = (InvTab4 *)((((ptrdiff_t)(xinvtab))+63L)&(-64L));

  /* Inversity of functions check:
     Check that myLinearSiteIndex3D is in fact the inverse
     of mySiteCoords3D, and that QDP_getSiteCoords is the
     inverse of QDP_linearSiteIndex()
  */
  for(p=0; p < 2; p++) {
    for(site=0; site < subgrid_vol_cb; site++) { 
      
      /* Linear site index */
      my_index = site + subgrid_vol_cb*p;
      QDP_getSiteCoords(gcoord, my_node, my_index);
      linear=QDP_getLinearSiteIndex(gcoord);

      if( linear != my_index ) { 
	printf("P%d cb=%d site=%d : QDP_getSiteCoords not inverse of QDP_getLinearSiteIndex(): my_index=%d linear=%d\n", my_node, p,site, my_index,linear);
      }

      mySiteCoords4D(gcoord, my_node, my_index);
      linear=myLinearSiteIndex4D(gcoord);

      if( linear != my_index ) { 
	printf("P%d cb=%d site=%d : mySiteCoords3D not inverse of myLinearSiteIndex3D(): my_index=%d linear=%d\n", my_node, p,site, my_index,linear);
      }
    }
  }


  /* Loop through sites - you can choose your path below */
  /* This is a checkerboarded order which is identical hopefully
     to QDP++'s rb2 subset when QDP++ is in a CB2 layout */
  for(p=0; p < 2; p++) { 
    for(t=0; t < subgrid_size[3]; t++) { 	   
      for(z=0; z < subgrid_size[2]; z++) {
	for(y=0; y < subgrid_size[1]; y++) { 
	  for(x=0; x < subgrid_size[0]/2; x++) {
	    
	    coord[0] = 2*x + p;
	    coord[1] = y;
	    coord[2] = z; 
	    coord[3] = t;
	  
	    /* Make global */
	    for(i=0; i < 4; i++) { 
	      coord[i] += subgrid_size[i]*node_coord[i];
	    }
	    
	    /* Index of coordinate -- NB this is not lexicographic
	       but takes into account checkerboarding in QDP++ */
	    qdp_index = QDP_getLinearSiteIndex(coord);
	    
	    /* Index of coordinate in my layout. -- NB this is not lexicographic
	       but takes into account my 3D checkerbaording */
	    my_index = myLinearSiteIndex4D(coord);
	    site_table[my_index] = qdp_index;

	    cb=parity(coord);
	    linear = my_index%subgrid_vol_cb;
	    
	    invtab[qdp_index].cb=cb;
	    invtab[qdp_index].linearcb=linear;
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
    for(site=0; site < subgrid_vol_cb; site++) {
      
      /* My local index */
      my_index = site + subgrid_vol_cb*p;
	
      /* Convert to QDP index */
      qdp_index = site_table[ my_index ];
      
      /* Switch QDP index to coordinates */
      QDP_getSiteCoords(gcoord, my_node,qdp_index);

      /* Convert back to cb3d index */
      linear = myLinearSiteIndex4D(gcoord);

      /* Check new cb,cbsite index matches the old cb index */
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
    for(site=0; site < subgrid_vol_cb; site++) {
      
      /* My local index */
      my_index = site + subgrid_vol_cb*p;
      mySiteCoords4D(gcoord, my_node, my_index);

      qdp_index = site_table[ my_index ];
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
  
  */
 
  /* This 4 is for the 4 tables: Table 1-4*/
  if ((shift_table = (int **)malloc(4*sizeof(int*))) == 0 ) {
    QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
    QMP_abort(1);
    
  }
  
  for(i=0; i < 4; i++) { 
    /* This 4 is for the 4 comms dierctions: */
    if ((shift_table[i] = (int *)malloc(4*subgrid_vol*sizeof(int))) == 0) {
      QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
      QMP_abort(1);
    }
  }

 
  /* Initialize the boundary counters */
  for(cb=0; cb < 2; cb++) {
    for(dir=0; dir < 4; dir++) {
      bound[cb][0][dir] = 0;	
      bound[cb][1][dir] = 0;	
      bound[cb][2][dir] = 0;	
      bound[cb][3][dir] = 0;	
    }
  }


  for(cb=0; cb < 2; cb++) { 
    for(site=0; site < subgrid_vol_cb; ++site) {
      
      index = cb*subgrid_vol_cb + site;
      
      /* Fetch site from site table */
      qdp_index = site_table[index];
      
      /* Get its coords */
      QDP_getSiteCoords(coord, my_node, qdp_index);
      
      /* Loop over directions building up shift tables */
      for(dir=0; dir < 4; dir++) {
	
	int fcoord[4], bcoord[4];
	int fnode, bnode;
	int blinear, flinear;
	
	/* Backwards displacement*/
	offs(bcoord, coord, dir, -1);
	bnode   = QDP_getNodeNumber(bcoord);
	blinear = QDP_getLinearSiteIndex(bcoord);

	/* Forward displacement */
	offs(fcoord, coord, dir, +1);
	fnode   = QDP_getNodeNumber(fcoord);
	flinear = QDP_getLinearSiteIndex(fcoord);

	/* Scatter:  decomp_{plus,minus} */
	/* Operation: a^F(shift(x,type=0),dir) <- decomp(psi(x),dir) */ 
	/* Send backwards - also called a receive from forward */
	if (bnode != my_node) {      
	  /* Offnode */
	  /* Append to Tail 1, increase boundary count */
	  /* This is the correct code */
	  shift_table[DECOMP_SCATTER][dir+4*index] 
	    = subgrid_vol_cb + bound[1-cb][DECOMP_SCATTER][dir];
	  
	  bound[1-cb][DECOMP_SCATTER][dir]++;
	  
	}
	else {                                           
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup */
	  shift_table[DECOMP_SCATTER][dir+4*index] = 
	    invtab[blinear].linearcb;
	}
	
	
	/* Scatter:  decomp_hvv_{plus,minus} */
	/* Operation:  a^B(shift(x,type=1),dir) <- U^dag(x,dir)*decomp(psi(x),dir) */
	/* Send forwards - also called a receive from backward */
	if (fnode != my_node) {
	  /* Offnode */
	  /* Append to Tail 1, increase boundary count */
	  shift_table[DECOMP_HVV_SCATTER][dir+4*index]           
	    = subgrid_vol_cb + bound[1-cb][DECOMP_HVV_SCATTER][dir];
	  
	  bound[1-cb][DECOMP_HVV_SCATTER][dir]++;                  

	}
	else {
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup */
	  shift_table[DECOMP_HVV_SCATTER][dir+4*index]           /* Onnode */
	    = invtab[flinear].linearcb ;
	}
	
	
	/* Gather:  mvv_recons_{plus,minus} */
	/* Operation:  chi(x) <-  \sum_dir U(x,dir)*a^F(shift(x,type=2),dir) */
	/* Receive from forward */
	if (fnode != my_node) {
	  /* Offnode */
	  /* Append to Tail 2, increase boundary count */

	  shift_table[RECONS_MVV_GATHER][dir+4*index] =
	    2*subgrid_vol_cb + (bound[cb][RECONS_MVV_GATHER][dir]);
	  
	  bound[cb][RECONS_MVV_GATHER][dir]++;

	}
	else {
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup. Note this is a recons post shift,
	     so the linear coordinate to invert is mine rather than the neighbours */
	  shift_table[RECONS_MVV_GATHER][dir+4*index] =
	    invtab[qdp_index].linearcb ;
	}
      
	/* Gather:  recons_{plus,minus} */
	/* Operation:  chi(x) +=  \sum_dir recons(a^B(shift(x,type=3),dir),dir) */
	/* Receive from backward */
	if (bnode != my_node) {

	  shift_table[RECONS_GATHER][dir+4*index] = 
	    2*subgrid_vol_cb + bound[cb][RECONS_GATHER][dir];
	  
	  bound[cb][RECONS_GATHER][dir]++;

	}
	else {
	  /* On node. Note the linear part of its (cb3, linear) bit,
	     using a reverse lookup. Note this is a recons post shift,
	     so the linear coordinate to invert is mine rather than the neighbours */
	  
	  shift_table[RECONS_GATHER][dir+4*index] = 
	    invtab[qdp_index].linearcb ;
	}
      } 
    }
  }

  /* Sanity check - make sure the sending and receiving counters match */
  for(cb=0; cb < 2; cb++) {
    for(dir=0; dir < 4; dir++) {

      /* Sanity 1: Must have same number of boundary sites on each cb for 
	 a given operation */
      for(i = 0; i < 4; i++) { 
	if (bound[1-cb][i][dir] != bound[cb][i][dir]) {
	  
	  QMP_error("SSE Wilson dslash - make_shift_tables: type 0 diff. cb send/recv counts do not match: %d %d",
		    bound[1-cb][i][dir],bound[cb][i][dir]);
	  QMP_abort(1);
	}
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
     
     Perhaps the best way to condsider this is to consider a value
     of shift_table[type][dir/site] that lands in the body. The
     shift table merely gives me a site index. But the data needs
     to be different for each direction for that site index. Hence 
     we need to replicate the body, for each dir. The 3xsubgrid_vol_cb
     is just there to take care of the buffers.

     Or another way to think of it is that there is a 'body element' index
     specified by the shift table lookup, and that dir is just the slowest
     varying index.
       
  */

  /* 4 dims, 4 types, rest of the magic is to align the thingie */
  xoffset_table = (halfspinor_array **)malloc(4*4*subgrid_vol*sizeof(halfspinor_array*)+63L);
  if( xoffset_table == 0 ) {
    QMP_error("init_wnxtsu3dslash: could not initialize offset_table[i]");
    QMP_abort(1);
  }
  /* This is the bit what aligns straight from AMD Manual */
  offset_table = (halfspinor_array**)((((ptrdiff_t)(xoffset_table)) + 63L) & (-64L));

  /* Walk through the shift_table and remap the offsets into actual
     pointers */

  /* DECOMP_SCATTER */
  num=0;
  for(dir =0; dir < Nd; dir++) { 

    /* Loop through all the sites. Remap the offsets either to local 
       arrays or pointers */
    offsite_found=0;
    for(site=0; site < subgrid_vol; site++) { 
      offset = shift_table[DECOMP_SCATTER][dir+4*site];
      if( offset >= subgrid_vol_cb ) { 
	/* Found an offsite guy. It's address must be to the send back buffer */
	/* send to back index = recv from forward index = 0  */
	offsite_found++;
	offset_table[ dir + 4*(site + subgrid_vol*DECOMP_SCATTER) ] =
	  send_bufs[0][num]+(offset - subgrid_vol_cb);
      }
      else { 
	/* Guy is onsite: This is DECOMP_SCATTER so offset to chi1 */
	offset_table[ dir + 4*(site + subgrid_vol*DECOMP_SCATTER) ] =
	chi1+shift_table[DECOMP_SCATTER][dir+4*site]+subgrid_vol_cb*dir;
      }
    }

    if( offsite_found > 0 ) { 
      /* If we found an offsite guy, next direction has to 
	 go into the next dir part of the send bufs */
      num++; 
    }
  }

  /* DECOMP_HVV_SCATTER */
  /* Restart num-s */
  num=0;
  for(dir =0; dir <Nd; dir++) { 
    offsite_found=0;
    for(site=0; site < subgrid_vol; site++) { 
      offset = shift_table[DECOMP_HVV_SCATTER][dir+4*site];
      if( offset >= subgrid_vol_cb ) { 
	/* Found an offsite guy. It's address must be to the send forw buffer */
	/* send to forward / receive from backward index = 1 */
	offsite_found++;

	offset_table[ dir + 4*(site + subgrid_vol*DECOMP_HVV_SCATTER) ] =
	  send_bufs[1][num]+(offset - subgrid_vol_cb);
      }
      else { 
	/* Guy is onsite. This is DECOMP_HVV_SCATTER so offset to chi2 */
	offset_table[ dir + 4*(site + subgrid_vol*DECOMP_HVV_SCATTER) ] =
	  chi2+shift_table[DECOMP_HVV_SCATTER][dir+4*site ]+subgrid_vol_cb*dir;
      }
    }
    if( offsite_found > 0 ) { 
      num++; 
    }
  }

  /* RECONS_MVV_GATHER */
  num=0;
  for(dir =0; dir <Nd; dir++) { 
    offsite_found=0;
    for(site=0; site < subgrid_vol; site++) { 
      offset = shift_table[RECONS_MVV_GATHER][dir+4*site];
      if( offset >= 2*subgrid_vol_cb ) { 
	/* Found an offsite guy. It's address must be to the recv from front buffer */
	/* recv_from front index = send to back index = 0 */
	offsite_found++;
	offset_table[ dir + 4*(site + subgrid_vol*RECONS_MVV_GATHER) ] =
	  recv_bufs[0][num]+(offset - 2*subgrid_vol_cb);
      }
      else { 
	/* Guy is onsite */
	/* This is RECONS_MVV_GATHER so offset with respect to chi1 */
	offset_table[ dir + 4*(site + subgrid_vol*RECONS_MVV_GATHER) ] =
	  chi1+shift_table[RECONS_MVV_GATHER][dir+4*site ]+subgrid_vol_cb*dir;
      }
    }
    if( offsite_found > 0 ) { 
      num++; 
    }
  }

  /* RECONS_GATHER */
  num=0;
  for(dir =0; dir <Nd; dir++) { 
    offsite_found=0;
    for(site=0; site < subgrid_vol; site++) { 
      offset = shift_table[RECONS_GATHER][dir+4*site];
      if( offset >= 2*subgrid_vol_cb ) { 
	/* Found an offsite guy. It's address must be to the recv from back buffer */
	/* receive from back = send to forward index =  1*/
	offsite_found++;
	offset_table[ dir + 4*(site + subgrid_vol*RECONS_GATHER) ] =
	  recv_bufs[1][num]+(offset - 2*subgrid_vol_cb);
      }
      else { 
	/* Guy is onsite */
	/* This is RECONS_GATHER so offset with respect to chi2 */
	offset_table[ dir + 4*(site + subgrid_vol*RECONS_GATHER ) ] = 
	chi2+shift_table[RECONS_GATHER][dir+4*site ]+subgrid_vol_cb*dir;
      }
    }
    if( offsite_found > 0 ) { 
      num++; 
    }
  }



  /* Free shift table - it is no longer needed. We deal solely with offsets */
  for(i=0; i < 4; i++) { 
    free( (shift_table)[i] );
  }
  free( shift_table );

  free( xinvtab );
  
  } 


void free_shift_tables(void) 
{
  free( xoffset_table );
  free( xsite_table );
}


#ifdef __cplusplus
}
#endif
