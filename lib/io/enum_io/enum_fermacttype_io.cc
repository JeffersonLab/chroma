#include "enum_fermacttype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace FermActTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theFermActTypeMap::Instance().registerPair(string("WILSON"), 
							   FERM_ACT_WILSON );

      success &=theFermActTypeMap::Instance().registerPair(string("UNPRECONDITIONED_WILSON"), 
					     FERM_ACT_UNPRECONDITIONED_WILSON);
      
	success &=theFermActTypeMap::Instance().registerPair(string("PARITY_BREAKING_WILSON"),
							     FERM_ACT_PARITY_BREAKING_WILSON);

	success &=theFermActTypeMap::Instance().registerPair(string("CLOVER"),
							     FERM_ACT_CLOVER);
	success &=theFermActTypeMap::Instance().registerPair(string("UNPRECONDITIONED_CLOVER"),FERM_ACT_UNPRECONDITIONED_CLOVER);

	success &=theFermActTypeMap::Instance().registerPair(string("DWF"),
							     FERM_ACT_DWF);

	success &=theFermActTypeMap::Instance().registerPair(string("UNPRECONDITIONED_DWF"),
							     FERM_ACT_UNPRECONDITIONED_DWF);

	success &=theFermActTypeMap::Instance().registerPair(string("PROJECTED_DWF"),
							     FERM_ACT_PROJECTED_DWF);

	success &=theFermActTypeMap::Instance().registerPair(string("OVERLAP_PARTIAL_FRACTION_4D"),
							     FERM_ACT_OVERLAP_PARTFRAC_4D );

	success &=theFermActTypeMap::Instance().registerPair(string("UNPRECONDITIONED_OVERLAP_CONTINUED_FRACTION_5D"),
							     FERM_ACT_UNPRECONDITIONED_OVERLAP_CONTFRAC_5D );

	success &=theFermActTypeMap::Instance().registerPair(string("ZOLOTAREV_5D"),
							     FERM_ACT_OVLAP_CONTFRAC_5D);

	success &=theFermActTypeMap::Instance().registerPair(string("OVERLAP_DWF"),
							     FERM_ACT_OVERLAP_DWF);

	success &=theFermActTypeMap::Instance().registerPair(string("EXTENDED_OVERLAP"),
							     FERM_ACT_EXTENDED_OVERLAP);

	success &=theFermActTypeMap::Instance().registerPair(string("UNPRECONDITIONED_EXTENDED_OVERLAP"), FERM_ACT_UNPRECONDITIONED_EXTENDED_OVERLAP);

	success &=theFermActTypeMap::Instance().registerPair(string("SMEARED_LAPLACIAN_WILSON"),
							     FERM_ACT_SMEARED_LAPLACIAN_WILSON);

	success &=theFermActTypeMap::Instance().registerPair(string("PLANAR_WILSON"),
							     FERM_ACT_PLANAR_WILSON);

	success &=theFermActTypeMap::Instance().registerPair(string("HAMBER_WU"),
							     FERM_ACT_HAMBER_WU);

	success &=theFermActTypeMap::Instance().registerPair(string("STAGGERED"),
							     FERM_ACT_STAGGERED);

	success &=theFermActTypeMap::Instance().registerPair(string("NAIK"),
							     FERM_ACT_NAIK);

	success &=theFermActTypeMap::Instance().registerPair(string("ASQTAD"),
							     FERM_ACT_ASQTAD);
      
      return success;
    }
    
    const string typeIDString = "FermActType";
    bool registered = registerAll();
  };

  using namespace FermActTypeEnv;

  //! Read an fermact type enum
  void read(XMLReader& xml_in,  const string& path, FermActType& t) {
    theFermActTypeMap::Instance().read(typeIDString,xml_in, path,t);
  }
  
  //! Write an fermact type  enum
  void write(XMLWriter& xml_out, const string& path, const FermActType& t) {
    theFermActTypeMap::Instance().write(typeIDString,xml_out, path, t);
  }
};
