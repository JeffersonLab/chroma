#include <shift_table_scalar.h>

#include <iostream>
#include <cstdlib>
#include <cstddef>
#include "allocate.h"
namespace CPlusPlusWilsonDslash { 

  namespace { 

    
  }



  ShiftTable::ShiftTable(const int latt_size[],  
			 void (*getSiteCoords)(int coord[], int node, int linearsite),
			 int (*getLinearSiteIndex)(const int coord[]),
			 int (*nodeNum)(const int coord[])
			 ) : Nd(4)
  {

    for(int mu=0; mu < 4; mu++)  {
      if ( latt_size[mu] % 2 != 0 ) {
	std::cerr << "This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even,  Your lattice is not like this: latt_size[" <<
	  mu <<"]="<<latt_size[mu]<< std::endl; 
	exit(1);
      }

      /* Copy this internally */
      tot_size[mu] = latt_size[mu];
    }

    /* Compute total volume and total checkerboarded volume */
    total_vol = tot_size[0];
    for(int mu=1; mu < Nd; mu++) { 
      total_vol *= tot_size[mu];
    }
    total_vol_cb = total_vol/2;


    //    int* inv_table = (int *)CPlusPlusWilsonDslash::alloc(total_vol*sizeof(int));
    // if( inv_table == (int *)NULL ) { 
    //   std::cerr << "Could not allocate site table " << std::endl;
    //   exit(1);
    // }
       

    /* Set up the path */
    /* This sets upt the inv table */
    // setupPathTable(getLinearSiteIndex, inv_table);

      
    xshift_table = (int *)CPlusPlusWilsonDslash::alloc(4*total_vol*2*sizeof(int)+Cache::CacheLineSize);
      
    if ( xshift_table == 0x0 ) {
      std::cerr << "Could not allocate xshift table" << std::endl;
      exit(1);
    }

    ptrdiff_t pad = 0;
    if ( (ptrdiff_t)xshift_table % Cache::CacheLineSize != 0 ) {
      pad=Cache::CacheLineSize-((ptrdiff_t)xshift_table % Cache::CacheLineSize);
    }
    shift_table = (int *)((char *)xshift_table + pad);


    xsite_table = (int *)CPlusPlusWilsonDslash::alloc(total_vol*sizeof(int)+Cache::CacheLineSize);
    
    if ( xsite_table == 0x0 ) {
      std::cerr << "Could not allocate site table " << std::endl;
      exit(1);
    }


    pad = 0;
    if ( (ptrdiff_t)xsite_table % Cache::CacheLineSize != 0 ) {
      pad=(ptrdiff_t)Cache::CacheLineSize-((ptrdiff_t)xsite_table % Cache::CacheLineSize);
    }
    site_table = (int *)((char *)xsite_table + pad);    

#pragma omp parallel for collapse(5)
    for(int p=0; p < 2; p++) { 	    
      for(int t=0; t < tot_size[3]; t++) { 
	for(int z=0; z < tot_size[2]; z++) {
	  for(int y=0; y < tot_size[1]; y++) {     
	    for(int x=0; x < tot_size[0]/2; x++) {
	      int coord[4];
	      
	      coord[0] = 2*x+p;
	      coord[1] = y;
	      coord[2] = z; 
	      coord[3] = t;
	      
	      /* Get the site and N-parity of the chosen victim */
	      int qdp_index = getLinearSiteIndex(coord); /* get the lexico index */
	      int my_index = myLinearSiteIndex4D(coord);
	      
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
#pragma omp parallel for collapse(2)
    for(int cb=0; cb < 2; cb++) {
      for(int site = 0; site < total_vol_cb; ++site) { 
        int fcoord[4], bcoord[4], coord[4];
	int blinear, flinear;
	
	int my_index = cb*total_vol_cb + site;
	
	int qdp_index = site_table[ my_index ];

	getSiteCoords(coord, 0, qdp_index); 
	
	for(int dir=0; dir < 4; dir++) {
	  
	  /* Backwards displacement*/
	  offs(bcoord, coord, dir, -1);
	  blinear = getLinearSiteIndex(bcoord);
	  
	  /* Forward displacement */
	  offs(fcoord, coord, dir, +1);
	  flinear = getLinearSiteIndex(fcoord);
	  
	  
	  /* Gather */
	  shift_table[dir+Nd*my_index ] = blinear;
	  shift_table[dir+Nd*(my_index+total_vol)] = flinear;
	}
      }
    }

    //    CPlusPlusWilsonDslash::dealloc(inv_table);
  }



  void ShiftTable::mySiteCoords4D(int gcoords[], int node, int linearsite)
  {
    int mu;
    int subgrid_cb_nrow[4];
    int tmp_coord[4];
    int cb,cbb;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = tot_size[mu];
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

  int ShiftTable::myLinearSiteIndex4D(const int gcoords[]) 
  {
    int mu;
    int subgrid_cb_nrow[4];
    int subgrid_cb_coord[4];
    int cb;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = tot_size[mu];
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

    return localSite4d(subgrid_cb_coord, subgrid_cb_nrow) + cb*total_vol_cb;
  }




#if 0
  // This subroutine allows us to block our lattice.
  // However, nonstandard blockings would require remapping the 
  // Input and output Fermions.
  // ************* THIS ROUTINE IS DISABLED *********************
  void ShiftTable::setupPathTable(int (*getLinearSiteIndex)(const int coord[]),
				  int* inv_table) 
    {
      // Blocking strategy: Block size in each dimension;
      int blockMaxSize[4] = {1,1,1,1};
      
      int nBlocksPerDim[4];
      int blockSizePerDim[4];
      
      
      for(int mu=0; mu < 4; mu++) {
	// Integer divide to get number of blocks in each dim
	nBlocksPerDim[mu] = tot_size[mu]/blockMaxSize[mu];
	if ( tot_size[mu] % blockMaxSize[mu] != 0 ) { 
	  nBlocksPerDim[mu]++;
	}
      }
      
      for(int mu=0; mu < 4; mu++) { 
	blockSizePerDim[mu] = blockMaxSize[mu] + tot_size[mu]%blockMaxSize[mu];
      }
      
      // Count all the blocks...
      int nBlocks = nBlocksPerDim[0];
      for(int mu=1; mu < 4; mu++){ 
	nBlocks *= nBlocksPerDim[mu];
      }
      
#if 0      
      std::cout << "There are " << nBlocks << " blocks" << std::endl;
      std::cout << "(Bx, By, Bz, Bt) = (" << nBlocksPerDim[0] 
	   << " , " << nBlocksPerDim[1]
	   << " , " << nBlocksPerDim[2]
	   << " , " << nBlocksPerDim[3] << ")" << std::endl;
#endif
      
      typedef struct { 
      int origin[4];
	int size[4];
      } Block;
      
      Block* blocklist = new Block[ nBlocks ];
      int blockpos=0;
      for(int t=0; t < nBlocksPerDim[3]; t++) {
	for(int z=0; z < nBlocksPerDim[2];z++) { 
	  for(int y=0; y < nBlocksPerDim[1]; y++){ 
	    for(int x=0; x < nBlocksPerDim[0]; x++) { 
	      
	      blocklist[ blockpos ].origin[0] = x*blockMaxSize[0];
	      blocklist[ blockpos ].origin[1] = y*blockMaxSize[1];
	      blocklist[ blockpos ].origin[2] = z*blockMaxSize[2];
	      blocklist[ blockpos ].origin[3] = t*blockMaxSize[3];
	      
	      blocklist[ blockpos ].size[0] = blockMaxSize[0];
	      blocklist[ blockpos ].size[1] = blockMaxSize[1];
	      blocklist[ blockpos ].size[2] = blockMaxSize[2];
	      blocklist[ blockpos ].size[3] = blockMaxSize[3];
	      
	      // Last x block
	      if( x == (nBlocksPerDim[0] - 1) ) {
		blocklist[ blockpos ].size[0] = tot_size[0]
		  -blocklist[blockpos].origin[0];
	      }
	      // Last y block
	      if( y == (nBlocksPerDim[1] - 1) ) {
		blocklist[ blockpos ].size[1] = tot_size[1]
		  -blocklist[blockpos].origin[1];
	      }
	      // Last z block
	      if( z == (nBlocksPerDim[2] - 1) ) {
		blocklist[ blockpos ].size[2] = tot_size[2]
		  -blocklist[blockpos].origin[2];
	      }
	      // Last t block
	      if( t == (nBlocksPerDim[3] - 1) ) {
		blocklist[ blockpos ].size[3] = tot_size[3]
		  -blocklist[blockpos].origin[3];
	      }
	      
	      blockpos++;
	    }
	  }
	}
      }
      
      
#if 0
      for(int i=0; i < nBlocks; i++){ 
	cout << "Block " << i << " : origin = ( "
	     << blocklist[i].origin[0] << " , "
	     << blocklist[i].origin[1] << " , "
	     << blocklist[i].origin[2] << " , " 
	     << blocklist[i].origin[3] << " ) \t size = ( "
	     << blocklist[i].size[0] << " , "
	     << blocklist[i].size[1] << " , "
	     << blocklist[i].size[2] << " , "
	     << blocklist[i].size[3] << " ) "  << std::endl;
      }
#endif
      
      /* Consistency check: If I add up all the sites in all the blocks
	 I should get the total number of sites on the lattice */
      int blocksites = 0;
      for(int i=0; i < nBlocks; i++) { 
	int n = blocklist[ i ].size[0]
	  * blocklist[i].size[1]
	  * blocklist[i].size[2]
	  * blocklist[i].size[3];
	blocksites += n;
      }
      
      if( blocksites != total_vol ) { 
	cerr << "Failed test: No of sites in all blocks = " << blocksites 
	     << " No of sites in lattice = " << total_vol << std::endl;
      }
      
      // Create paths -- do a lexicographic walk in each block
      // and bin the points into even-odd path arrays...
      
      // First I have to make the even-odd path arrays...
      xpath_table =(int *)CPlusPlusWilsonDslash::alloc(sizeof(int)*total_vol+Cache::CacheLineSize );
      if( xpath_table == (int *)NULL ) { 
	cerr << "Failed to allocate xpath_table" << std::endl;
	exit(1);
      }
      ptrdiff_t pad = 0;
      if ( (ptrdiff_t)xpath_table % Cache::CacheLineSize != 0 ) {
	pad=(ptrdiff_t)Cache::CacheLineSize-((ptrdiff_t)xpath_table % Cache::CacheLineSize);
      }
      path_table = (int *)((char *)xpath_table + pad);
      
      // Go through blocklist
      int index[2]; index[0]=0; index[1]=0;
      
      for(int i=0; i < nBlocks; i++) { 
	for(int t=0; t < blocklist[i].size[3]; t++) { 
	  for(int z=0; z < blocklist[i].size[2]; z++) { 
	    for(int y=0; y < blocklist[i].size[1]; y++) { 
	      for(int x=0; x < blocklist[i].size[0]; x++) { 
		
		int coord[4];
		coord[0] = blocklist[i].origin[0]+x;
		coord[1] = blocklist[i].origin[1]+y;
		coord[2] = blocklist[i].origin[2]+z;
		coord[3] = blocklist[i].origin[3]+t;
		
		int cb = parity(coord);
		int my_index = cb*total_vol_cb + index[cb];

		// path_table[my_index] works with QDP++ indices
		path_table[ my_index ] = getLinearSiteIndex(coord);
		index[cb]++;

		// inv_table takes QDP indices to my index scheme
		inv_table[ path_table[my_index] ] = my_index;
	      }
	    }
	  }
	}
      }
      
      delete [] blocklist;
    }
#endif

}
