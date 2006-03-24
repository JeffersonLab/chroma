// $Id: szin_read_obj_funcmap.cc,v 2.3 2006-03-24 22:26:32 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/szin_read_obj_funcmap.h"
#include "io/szin_io.h"
#include "io/readszin.h"
#include "io/readszinqprop_w.h"

namespace Chroma
{
 
  //! IO function map environment
  /*! \ingroup inlineio */
  namespace SZINReadObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //! Read a propagator
      void SZINReadLatProp(const string& buffer_id,
			   const string& file)
      {
	LatticePropagator obj;
	XMLReader record_xml;
    
	readSzinQprop(record_xml, obj, file);

	XMLBufferWriter file_xml;
	push(file_xml, "SZIN");
	pop(file_xml);

	TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }

      //! Read a gauge field in floating precision
      void SZINReadArrayLatColMat(const string& buffer_id,
				  const string& file)
      {
	XMLReader record_xml;
	multi1d<LatticeColorMatrix> obj(Nd);

	readSzin(record_xml, obj, file);

	XMLBufferWriter file_xml;
	push(file_xml, "SZIN");
	pop(file_xml);

	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(buffer_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }

    }  // end anonymous namespace


    bool registerAll(void) 
    {
      bool success = true;
      success &= TheSZINReadObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								    SZINReadLatProp);
      success &= TheSZINReadObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								    SZINReadArrayLatColMat);
      return success;
    }

    bool registered = registerAll();
  }

}
