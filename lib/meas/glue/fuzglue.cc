/*! \file
 *  \brief Compute 'fuzzy' (blocked) glueball correlation functions
 */

#include "chromabase.h"
#include "meas/glue/block.h"
#include "meas/glue/fuzglue.h"
#include "meas/glue/gluecor.h"
#include "meas/glue/polycor.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{ 
  //! Compute 'fuzzy' (blocked) glueball correlation functions
  /*! 
   * \ingroup glue
   *
   * Driver for computation of 'fuzzy' (blocked) glueball correlation functions
   * using Tepers 'fuzzying' method.
   *
   * Warning: this works only for Nd = 4 ! (Construction of glueball states)
   * Warning: this works only for Nc = 2 and 3 ! (Projection of blocked links)
   *
   * \param xml_out    xml file object ( Write )
   * \param xml_group  string used for writing xml data ( Read )
   * \param u          gauge field ( Read )
   * \param j_decay    direction along which the exponentrial is allowed to decay ( Read )
   * \param BlkAccu    accuracy in fuzzy link projection ( Read )
   * \param BlkMax     maximum number of iterations in fuzzy link projection ( Read ) 
   */

  void fuzglue(XMLWriter& xml_out, const string& xml_group,
	       const multi1d<LatticeColorMatrix>& u, 
	       int j_decay, const Real& BlkAccu, int BlkMax)
  {
    START_CODE();

    multi1d<LatticeColorMatrix> u_fuz(Nd);
    multi1d<LatticeColorMatrix> u_tmp(Nd);

    int bl_level_max;
    int bl_level;
    int block_latt;
    int mu;

    if( Nd != 4 )
    {
      QDPIO::cout << " Glueball construction works only for Nd = 4 !" << endl;
      END_CODE();
      return;
    }

    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    /* Determine the maximum blocklevel possible orthogonal to direction j_decay */
    bl_level_max = QDP::Layout::lattSize()[0];
    for(mu = 0; mu < Nd; ++mu)
    {
      if ( mu != j_decay )
      {
	block_latt = QDP::Layout::lattSize()[mu];
	bl_level = 0;
	while( (block_latt > 2) && ((block_latt & 1) == 0) )
	{
	  block_latt = block_latt / 2;
	  bl_level = bl_level + 1;
	}
	if ( bl_level < bl_level_max )
	  bl_level_max = bl_level;
      }
    }

  
    /* Copy u's to u_fuz */
    u_fuz = u;

    /* Do unblocked measurements */
    bl_level = 0;
    gluecor(xml_out, "GlueCorr", u_fuz, phases, bl_level);
    polycor(xml_out, "PolyCorr", u_fuz, phases, bl_level);

    /* Loop over blocking levels */
    while( bl_level < bl_level_max )
    {
      /* Do the blocking of links perpendicular to direction j_decay */
      for(mu = 0; mu < Nd; ++mu)
      {
	if( mu != j_decay )
	  block(u_tmp[mu], u_fuz, mu, bl_level, BlkAccu, BlkMax, j_decay);
      }

      u_fuz = u_tmp;
    
      /* Advance bl_level, after blocking is complete */
      bl_level = bl_level + 1;

      /* Do blocked measurements */
      gluecor(xml_out, "GlueCorr", u_fuz, phases, bl_level);
      polycor(xml_out, "PolyCorr", u_fuz, phases, bl_level);
    }

    END_CODE();
  }

}
