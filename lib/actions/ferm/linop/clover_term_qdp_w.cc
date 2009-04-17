// $Id: clover_term_qdp_w.cc,v 3.22 2009-04-17 02:05:32 bjoo Exp $
/*! \file
 *  \brief Clover term linear operator
 *
 *  This particular implementation is specific to a scalar-like
 *  architecture
 */

#include "chromabase.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "meas/glue/mesfield.h"


namespace Chroma 
{ 
#if 0
  // Reader/writers
  void read(XMLReader& xml, const string& path, PrimitiveClovTriang& param)
  {
    QDPIO::cerr << __func__ << ": clover term reader not implemented" << endl;
    QDP_abort(1);
  }

  void write(XMLWriter& xml, const string& path, const PrimitiveClovTriang& d)
  {
    xml.openTag(path);

    XMLWriterAPI::AttributeList alist;

    xml.openTag("Diag");
    for(int i=0; i < 2; ++i)
    {
      for(int j=0; j < 2*Nc; ++j)
      {
	alist.clear();
	alist.push_back(XMLWriterAPI::Attribute("block", i));
	alist.push_back(XMLWriterAPI::Attribute("col", j));

	xml.openTag("elem", alist);
	xml << d.diag[i][j];
	xml.closeTag();
      }
    }
    xml.closeTag(); // Diag

    xml.openTag("Offd");
    for(int i=0; i < 2; ++i)
    {
      for(int j=0; j < 2*Nc*Nc-Nc; ++j)
      {
	alist.clear();
	alist.push_back(XMLWriterAPI::Attribute("block", i));
	alist.push_back(XMLWriterAPI::Attribute("col", j));

	xml.openTag("elem", alist);
	xml << d.offd[i][j];
	xml.closeTag();
      }
    }
    xml.closeTag(); // Offd
    xml.closeTag(); // path
  }
#endif



} // Namespace Chroma



