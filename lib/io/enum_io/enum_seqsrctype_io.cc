#include "enum_seqsrctype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace SeqSourceTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theSeqSourceTypeMap::Instance().registerPair(string("NUCL_U_UNPOL"), 
							     SEQ_SOURCE_TYPE_NUCL_U_UNPOL );

      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_D_UNPOL"), 
							     SEQ_SOURCE_TYPE_NUCL_D_UNPOL  );
      
      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_U_POL"),
							     SEQ_SOURCE_TYPE_NUCL_U_POL );

      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_D_POL"),
							       SEQ_SOURCE_TYPE_NUCL_D_POL);
      
      success &=theSeqSourceTypeMap::Instance().registerPair(string("DELTA_U_UNPOL"),
							     SEQ_SOURCE_TYPE_DELTA_U_UNPOL );
      
      success &=theSeqSourceTypeMap::Instance().registerPair(string("DELTA_D_UNPOL"),
							     SEQ_SOURCE_TYPE_DELTA_D_UNPOL );
      
      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_U_UNPOL_NONREL"),
							       SEQ_SOURCE_TYPE_NUCL_U_UNPOL_NONREL );
      
      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_D_UNPOL_NONREL"),
							     SEQ_SOURCE_TYPE_NUCL_D_UNPOL_NONREL);

      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_U_POL_NONREL"),
							     SEQ_SOURCE_TYPE_NUCL_U_POL_NONREL);

      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_D_POL_NONREL"),
							     SEQ_SOURCE_TYPE_NUCL_D_POL_NONREL);
      
      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_U_MIXED_NONREL"),
							     SEQ_SOURCE_TYPE_NUCL_U_MIXED_NONREL);

      success &=theSeqSourceTypeMap::Instance().registerPair(string("NUCL_D_MIXED_NONREL"),   
							     SEQ_SOURCE_TYPE_NUCL_D_MIXED_NONREL);

      success &=theSeqSourceTypeMap::Instance().registerPair(string("PION"),
							       SEQ_SOURCE_TYPE_PION);
      
      return success;
    }

    const string typeIDString = "SeqSourceType";
    const bool registered = registerAll();
  };
  using namespace SeqSourceTypeEnv;

  //! Read a sequential source type enum
  void read(XMLReader& xml_in,  const string& path, SeqSourceType& t) {
    theSeqSourceTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write a sequential source  type enum
  void write(XMLWriter& xml_out, const string& path, const SeqSourceType& t) {
    theSeqSourceTypeMap::Instance().write(typeIDString,xml_out, path, t);
  }
};
