// $Id: seqsrc_funcmap_w.cc,v 1.3 2005-03-18 05:12:37 edwards Exp $
/*! \file
 *  \brief Sequential source function map
 */

#include "meas/hadron/seqsrc_funcmap_w.h"
#include "meas/hadron/barseqsrc_w.h"
#include "meas/hadron/mesonseqsrc_w.h"

namespace Chroma
{
 
  namespace SeqSourceCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_UNPOL"), 
								  barNuclUUnpol);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_UNPOL"), 
								  barNuclDUnpol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_POL"),
								  barNuclUPol);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_POL"),
								  barNuclDPol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("DELTA_U_UNPOL"),
								  barDeltaUUnpol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("DELTA_D_UNPOL"),
								  barDeltaDUnpol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_UNPOL_NONREL"),
								  barNuclUUnpolNR);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_UNPOL_NONREL"),
								  barNuclDUnpolNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_POL_NONREL"),
								  barNuclUPolNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_POL_NONREL"),
								  barNuclDPolNR);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_MIXED_NONREL"),
								  barNuclUMixedNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_MIXED_NONREL"),   
								  barNuclDMixedNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A0_1"),
								  mesPionA01SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_X_1"),
								  mesPionRhoX1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Y_1"),
								  mesPionRhoY1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-B1_Z_1"),
								  mesPionB1Z1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Z_1"),
								  mesPionRhoZ1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-B1_Y_1"),
								  mesPionRhoZ1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-B1_X_1"),
								  mesPionRhoZ1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-PION_2"),
								  mesPion1Pion2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A0_2"),
								  mesPionA02SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_X_2"),
								  mesPionRhoX2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_X_2"),
								  mesPionRhoX2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Y_2"),
								  mesPionRhoY2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A1_Z_1"),
								  mesPionA1Z1SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Z_2"),
								  mesPionRhoZ2SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A1_Y_1"),
								  mesPionA1Y1SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A1_X_1"),
								  mesPionA1X1SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION"),
								  mesPion1Pion1SeqSrc);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION_1"), // same as PION
								  mesPion1Pion1SeqSrc);
      
      return success;
    }

    bool registered = registerAll();
  };

}
