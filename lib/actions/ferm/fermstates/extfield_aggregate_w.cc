// $Id: extfield_aggregate_w.cc,v 1.4 2007-10-30 02:56:02 kostas Exp $
/*! \file
 *  \brief External field aggregate
 */

#include "chromabase.h"

#include "actions/ferm/fermstates/extfield_factory_w.h"
#include "actions/ferm/fermstates/extfield_aggregate_w.h"

#include "state.h"
#include "create_state.h"


namespace Chroma 
{

  //! External fields
  /*! \ingroup fermstates */
  namespace ExternalFieldEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! Construct zero
      ExternalField* zeroFunc(XMLReader& xml_in,
			      const std::string& path)
      {
	return new ZeroExternalField();
      }

      //! Construct constant magnetic field
      ExternalField* ConstantMagneticFunc(XMLReader& xml_in,
					  const std::string& path)
      {
	ConstantMagneticParams p ;
	read(xml_in,path,p);
	return new ConstantMagneticExternalField(p);
      }

#if 0
      //! Construct linear term
      ExternalField* linearFunc(XMLReader& xml_in,
				const std::string& path)
      {
	return new LinearExternalField(LinearParams(xml_in, path));
      }
#endif
      
    } // end anonymous namespace

    

  //! Reader
  /*! @ingroup sources */
    void read(XMLReader& xml, const string& path, ExternalFieldEnv::ConstantMagneticParams& param){
      
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
	{
	case 1:
	  break;
	  
	default:
	  QDPIO::cerr << __func__ << ": parameter version " << version 
		      << " unsupported." << endl;
	  QDP_abort(1);
	}

      read(paramtop, "time_dir",  param.time_dir);
      read(paramtop, "Bfield", param.Bfield);
      read(paramtop, "dir", param.dir);

    }
  
  //! Writer
  /*! @ingroup sources */
    void write(XMLWriter& xml, const string& path, const ExternalFieldEnv::ConstantMagneticParams& param){
 
      push(xml,path);
      int version;
      write(xml, "version", version);
      write(xml, "time_dir",  param.time_dir);
      write(xml, "Bfield", param.Bfield);
      write(xml, "dir", param.dir);
      pop(xml);
    }


    multi1d< Handle< ExternalField > > reader(XMLReader& xml, 
					      const std::string& path){
      multi1d< Handle< ExternalField > > ext_field(Nd) ;
       XMLReader paramtop(xml, path);
       string name ;
       read(paramtop, "Name", name);
       /**
              Handle< ExternalField > ef = Chroma::TheExternalFieldFactory::Instance().createObject(paramtop,name) ;
      for(int mu(0);mu<Nd;mu++){
	ext_field[mu] = ef ;
      }
       **/
       QDPIO::cout<<"Found external field: "<<name<<endl ;

      return ext_field ;
    }

    // Construct linear function
    LatticeComplex
    ZeroExternalField::operator()(int dir) const
    {
      START_CODE();

      LatticeComplex d = 1;

      END_CODE();
      return d;
    }

    // Construct constant magnetic field. Needs fix up...
    LatticeComplex
    ConstantMagneticExternalField::operator()(int mu) const
    {
      START_CODE();

      LatticeComplex U(1.0) ;
      Real A(0.5*Bfield) ;
      Complex C = cmplx(cos(A),sin(A));
      if(Nd==4){
	multi1d<int> dd(2) ;
	for( int d(0);d<Nd;d++){
	  if((d!=time_dir)&&(d!=dir)&&(d!=mu)){
	    A *= epsilon(d,mu,dir) ;
	    LatticeReal B = A*Layout::latticeCoordinate(mu) ; 
	    U = cmplx(cos(B),sin(B));
	    return U ;
	  }
	}
      }
      else//do nothing
	QDPIO::cerr<<"OOOPS! don't know how to set up magnetic field in other than 4 dimenssions\n"; 
      
      END_CODE();
      return U;
    }



    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Chroma::TheExternalFieldFactory::Instance().registerObject(string("ZERO"), zeroFunc);
	QDPIO::cerr<<"registered ZERO external field\n";

	success &= Chroma::TheExternalFieldFactory::Instance().registerObject(string("CONSTANT_MAGNETIC"), ConstantMagneticFunc);
	QDPIO::cerr<<"registered CONSTANT_MAGNETIC field\n";

	registered = true;
      }
      return success;
    }
  }  // end namespace ExternalFieldEnv

}  // end namespace Chroma



  
