// $Id: seqsrc_funcmap_w.cc,v 1.1 2005-03-07 02:55:20 edwards Exp $
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

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION"),
								  mesPionSeqSrc);
      
      return success;
    }

    bool registered = registerAll();
  };

}
