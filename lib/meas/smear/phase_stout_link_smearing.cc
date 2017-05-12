/*! \file
 *  \brief Stout link smearing
 */

#include "chromabase.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/phase_stout_link_smearing.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, PhaseStoutLinkSmearingEnv::Params& param)
  {
    PhaseStoutLinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const std::string& path, const PhaseStoutLinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace PhaseStoutLinkSmearingEnv
  {
    namespace
    {
      //! Callback function
      LinkSmearing* createSource(XMLReader& xml_in,
				 const std::string& path)
      {
	return new LinkSmear(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name = "PHASE_STOUT_SMEAR";
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinkSmearingFactory::Instance().registerObject(name, createSource);
	registered = true;
      }
      return success;
    }


    //! Parameters for running code
    Params::Params(XMLReader& xml, const std::string& path) : StoutLinkSmearingEnv::Params(xml,path)
    {
      XMLReader paramtop(xml, path);
      read(paramtop, "k", k);
      read(paramtop, "zeta", zeta);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "version", version);
      write(xml, "LinkSmearingType", PhaseStoutLinkSmearingEnv::name);
      write(xml, "link_smear_num", link_smear_num);
      write(xml, "link_smear_fact", link_smear_fact);
      write(xml, "smear_dirs", smear_dirs);
      write(xml,"k",k);
      write(xml,"zeta",zeta);
      pop(xml);
    }


    //! Smear the links
    void
    LinkSmear::operator()(multi1d<LatticeColorMatrix>& u) const
    {
      // Now stout smear
      StoutLinkSmearingEnv::LinkSmear ss(params) ;
      ss(u);
      for(int d(0);d<Nd;d++)
	if(params.smear_dirs[d]){
	  if(params.k[d]!=0){
	    QDPIO::cout<<" Adding phase to direction: "<<d<<std::endl ;
	    Real f = twopi / Real(Layout::lattSize()[d])*
	      params.zeta*Real(params.k[d]) ;
	    Complex c = cmplx(cos(f),sin(f)) ;
	    QDPIO::cout<<"    exp(i*phase)= "<<c<<std::endl ;
	    QDPIO::cout<<"    2*pi= "<<twopi<<std::endl ;
	    u[d]*=c ;
	  }
	}
    }
  }
}

  
