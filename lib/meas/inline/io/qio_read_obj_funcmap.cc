// $Id: qio_read_obj_funcmap.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#include "named_obj.h"
#include "meas/inline/io/qio_read_obj_funcmap.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/eigeninfo.h"

namespace Chroma
{
 
  //! IO function map environment
  /*! \ingroup inlineio */
  namespace QIOReadObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //! Read a propagator
      void QIOReadLatProp(const string& buffer_id,
			  const string& file, 
			  QDP_serialparallel_t serpar)
      {
	LatticePropagator obj;
	XMLReader file_xml, record_xml;

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }


      //! Read a single prec propagator
      void QIOReadLatPropF(const string& buffer_id,
			   const string& file, 
			   QDP_serialparallel_t serpar)
      {
	LatticePropagatorF obj;
	XMLReader file_xml, record_xml;

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }


      //! Read a double prec propagator
      void QIOReadLatPropD(const string& buffer_id,
			   const string& file, 
			   QDP_serialparallel_t serpar)
      {
	LatticePropagatorD obj;
	XMLReader file_xml, record_xml;

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }



      //! Read a fermion
      void QIOReadLatFerm(const string& buffer_id,
			  const string& file, 
			  QDP_serialparallel_t serpar)
      {
	LatticeFermion obj;
	XMLReader file_xml, record_xml;

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create<LatticeFermion>(buffer_id);
	TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }


#if 0
      // RGE: FOR SOME REASON, QDP CANNOT CAST A DOUBLE TO FLOATING HERE. NEED TO FIX.

      //! Read a single prec fermion
      void QIOReadLatFermF(const string& buffer_id,
			   const string& file, 
			   QDP_serialparallel_t serpar)
      {
	LatticeFermionF obj;
	XMLReader file_xml, record_xml;

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create<LatticeFermion>(buffer_id);
	TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }


      //! Read a double prec fermion
      void QIOReadLatFermD(const string& buffer_id,
			   const string& file, 
			   QDP_serialparallel_t serpar)
      {
	LatticeFermionD obj;
	XMLReader file_xml, record_xml;

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create<LatticeFermion>(buffer_id);
	TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }
#endif

      //! Read a propagator
      void QIOReadArrayLatColMat(const string& buffer_id,
				 const string& file, 
				 QDP_serialparallel_t serpar)
      {
	multi1d<LatticeColorMatrix> obj;
	XMLReader file_xml, record_xml;

	obj.resize(Nd);    // BAD BAD BAD - FIX THIS

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(buffer_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }

      void QIOReadEigenInfo(const string& buffer_id,
			    const string& file,
			    QDP_serialparallel_t serpar)
      {

      
	multi1d<Real64> evalsD;

	// BUG? Need to read these as an array even though there is only one
	// and also I need to know in advance the array size.
	multi1d<Real64> largestD(1);

	// Would like this to be double, but casting double prec to single prec doesn't work.
	LatticeFermion evecD;

	// Stuff 
	XMLReader file_xml, record_xml, largest_record_xml;

	// Open file
	QDPFileReader to(file_xml,file,serpar);

	// Create object
	TheNamedObjMap::Instance().create< EigenInfo >(buffer_id);

	// Read largest ev plus XML file containing number of eval/evec pairs
	read(to, largest_record_xml, largestD);

	// Extract number of EVs from XML
	int size;
	try { 
	  read(largest_record_xml, "/NumElem/n_vec", size);
	}
	catch(const std::string& e) { 
	  QDPIO::cerr<< "Caught Exception while reading XML: " << e << endl;
	  QDP_abort(1);
	}

	// Set largest EV
	TheNamedObjMap::Instance().getData< EigenInfo >(buffer_id).getLargest()=largestD[0];

	// REsize eval arrays so that IO works correctly
	evalsD.resize(size);

	// Read evals
	read(to, record_xml, evalsD);

	QDPIO::cout << "Read " << evalsD.size() << "Evals " << endl;
	for(int i=0; i < evalsD.size(); i++) { 
	  QDPIO::cout << "Eval["<<i<<"] = " << Real(evalsD[i]) << endl;
	}

	// Resize eval array 
	TheNamedObjMap::Instance().getData< EigenInfo >(buffer_id).getEvalues().resize(evalsD.size());
      
	// Downcast evals to Real() precision
	for (int i=0; i<evalsD.size(); i++)
	  TheNamedObjMap::Instance().getData< EigenInfo >(buffer_id).getEvalues()[i] = Real(evalsD[i]);


	// Resize evec arrays
	TheNamedObjMap::Instance().getData< EigenInfo >(buffer_id).getEvectors().resize(evalsD.size());

	// Loop and read evecs
	for (int i=0; i<evalsD.size(); i++)
	{
	  XMLReader record_xml_dummy;
	  read(to, record_xml_dummy, evecD);
	  LatticeFermion evecR;
	  evecR=evecD;

	  (TheNamedObjMap::Instance().getData< EigenInfo >(buffer_id).getEvectors())[i]=evecR;
	}
 
	// Set File and Record XML throw away dummy XMLs
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);

	// Done - That too was unnecessarily painful
	close(to);
      }
    }  // end namespace


    bool registerAll(void) 
    {
      bool success = true;
      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								   QIOReadLatProp);
      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorF"), 
								   QIOReadLatPropF);
      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorD"), 
								   QIOReadLatPropD);

      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								   QIOReadLatFerm);
//      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermionF"), 
//								   QIOReadLatFermF);
//      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermionD"), 
//								   QIOReadLatFermD);

      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								   QIOReadArrayLatColMat);

      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("EigenInfo"),
								   QIOReadEigenInfo);

      return success;
    }

    bool registered = registerAll();
  }

}
