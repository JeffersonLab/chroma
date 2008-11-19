// -*- C++ -*-
// $Id: inline_stoch_hadron_w.h,v 1.8 2008-11-19 18:58:53 kostas Exp $
/*! \file
 * \brief Inline measurement of stochastic hadron operator (mesons and baryons).
 *
 * spectroscopy
 */

#ifndef __inline_stoch_hadron_h__
#define __inline_stoch_hadron_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
//#include <map>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochHadronEnv 
  {
    extern const std::string name;
    bool registerAll();
    
    // The flavors
    enum Flavor {up, down, strange, charm, bottom};

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
       
      unsigned long      frequency;
    
      struct Param_t
      {
	int          mom2_max;              /*!< (mom)^2 <= mom2_max */
	
	//Operators can be mesons or baryons (quark-antiquark or three quarks).
	multi1d<GroupXML_t> ops ; /*!< array with hadron operators xml */
	//Each operator needs its own output file
	//should be defined in the GroupXML_t
      
	GroupXML_t  smearing; /*!< xml holding smearing params. only one allowed */
	// for the moment ignore displacements...
	multi1d<GroupXML_t>  displace; /*!< xml holding displacement params */

	GroupXML_t   link_smear;         /*!< link smearing xml one for all*/

	multi1d<GroupXML_t> quarks ;     /*! dilutions */
	
      } param;
    
    
      struct NamedObject_t
      {
	std::string         gauge_id;
	//std::string         meson_db;
	//std::string         baryon_db;
      } named_obj;
      
      std::string xml_file;  // Alternate XML file pattern

      void write(XMLWriter& xml_out, const std::string& path);

    };

    struct MesonOp{
      int g ;
      string file;
    } ;

    struct BaryonOp{
      int g ;
      string file;
    } ;
  

    /**
    void meson(DComplex& corr,
	       const GroupXML_t& grpXML,
	       const LatticeComplex& phase,
	       const LatticeFermion& eta,
	       const LatticeFermion& chi,
	       const Subset& s) ; 
    
    void baryon(DComplex& corr,
		const GroupXML_t& grpXML,
		const LatticeComplex& phase,
		const LatticeFermion& eta1,
		const LatticeFermion& eta2,
		const LatticeFermion& eta3,
		const Subset& s) ;
    
    **/
    void ParseMeson(MesonOp&,const GroupXML_t&) ;
    void ParseBaryon(BaryonOp&,const GroupXML_t&) ;

  //! Inline measurement of stochastic baryon operators
  /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement{
    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 
      
    private:
 

      Params params;
      /**
      map<string, void (*)(DComplex& , 
			   const GroupXML_t& ,  
			   const LatticeComplex& , 
			   const LatticeFermion& , 
			   const LatticeFermion& ,  
			   const Subset& )> mesons ;

      map<string, void (*)(DComplex& ,
			   const GroupXML_t& ,
			   const LatticeComplex& ,
			   const LatticeFermion& ,
			   const LatticeFermion& ,
			   const LatticeFermion& ,
			   const Subset& )> baryons ;
      **/
      map<string, void (*)(MesonOp&,const GroupXML_t& )> mesons ;

      map<string, void (*)(BaryonOp&, const GroupXML_t&)> baryons ;

      void setUpMaps(){
	mesons["PION"]  = &ParseMeson  ;
	baryons["NUCLEON"] = &ParseBaryon ;
      }
      
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {setUpMaps();}
      InlineMeas(const InlineMeas& p) : params(p.params) {setUpMaps();}
      
      unsigned long getFrequency(void) const {return params.frequency;}
      
      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 
      
    };

  }; // name space InlineHadronEnv

};

#endif
