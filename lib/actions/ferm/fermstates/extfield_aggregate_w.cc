// $Id: extfield_aggregate_w.cc,v 1.11 2009-05-04 17:11:50 caubin Exp $
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

        //construct the antisymmetric tensor
    int epsilon(int i, int j, int k){
      if( (i==j)||(j==k)|| (k==i) )
	return 0 ;
      if( ((i<j)&&(j<k)) || ((j<k)&&(k<i)) || ((k<i)&&(i<j)) )
	return 1 ;
      else 
	return -1 ;
    }


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
      read(paramtop, "t_dir",  param.t_dir);
      read(paramtop, "y_dir",  param.y_dir);
      read(paramtop, "b_dir",  param.b_dir);
      read(paramtop, "Bfield", param.Bfield);
      if(paramtop.count("patch") !=0 )
	read(paramtop, "patch", param.patch);
      else
	param.patch = false ;
    }
  
  //! Writer
  /*! @ingroup sources */
    void write(XMLWriter& xml, const string& path, const ExternalFieldEnv::ConstantMagneticParams& param){
 
      push(xml,path);
      int version;
      write(xml, "version", version);
      write(xml, "t_dir",  param.t_dir);
      write(xml, "y_dir",  param.y_dir);
      write(xml, "b_dir",  param.b_dir);
      write(xml, "Bfield", param.Bfield);
      write(xml, "patch", param.patch);
      pop(xml);
    }


    //! Reader
    /*! @ingroup sources */
    void read(XMLReader& xml, const string& path, ExternalFieldEnv::LinearElectricParams& param){
      
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
      read(paramtop, "t_dir",  param.t_dir);
      read(paramtop, "x_dir",  param.x_dir);
      read(paramtop, "x_src",  param.x_src);
      read(paramtop, "Efield", param.Efield);
      if(paramtop.count("patch") !=0 )
        read(paramtop, "patch", param.patch);
      else
        param.patch = false ;
    }
  
  //! Writer
  /*! @ingroup sources */
    void write(XMLWriter& xml, const string& path, const ExternalFieldEnv::LinearElectricParams& param){
 
      push(xml,path);
      int version;
      write(xml, "version", version);
      write(xml, "t_dir",  param.t_dir);
      write(xml, "x_dir",  param.x_dir);
      write(xml, "x_src",  param.x_src);
      write(xml, "Efield", param.Efield);
      write(xml, "patch", param.patch);
      pop(xml);
    }


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
	ExternalFieldEnv::read(xml_in,path,p);
	return new ConstantMagneticExternalField(p);
      }

      
      //! Construct linear electric field
      ExternalField* LinearElectricFunc(XMLReader& xml_in,
					const std::string& path)
      {
	LinearElectricParams p ;
	ExternalFieldEnv::read(xml_in,path,p);
	return new LinearElectricExternalField(p);
      }
	  
      
    } // end anonymous namespace

    

    Handle< ExternalField > reader(XMLReader& xml, 
				   const std::string& path){
      XMLReader paramtop(xml, path);
      
      string name ;
      string ext_field_path = "ExternalField";
      
      read(paramtop, ext_field_path+"/Name", name);
      QDPIO::cout<<"Found external field: "<<name<<endl ;

      //Handle< ExternalField > ef ; // just to make the code compile

      /* THE FOLLOWING BIT DOES NOT WORK... why??? NEED FIX*/
      Handle< ExternalField > ef(Chroma::TheExternalFieldFactory::Instance().createObject(name,paramtop,ext_field_path)) ;
       /**/

      return ef ;
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

      LatticeComplex U ;
      if(mu!=t_dir){
	Real A;
	if(Nd==4){
	  if(mu==y_dir){
	    A = epsilon(x_dir,y_dir,b_dir)*Bfield ;
	    QDPIO::cout<<__func__<<" A_"<<mu<<"="<<A<<endl ;
	    LatticeReal B = A*Layout::latticeCoordinate(x_dir) ; 
	    U = cmplx(cos(B),sin(B));
	    return U ;
	  }
	  if((mu==x_dir)&&(patch)){// do the patch 
	    LatticeBoolean XX = (Layout::latticeCoordinate(x_dir)==
				(Layout::lattSize()[x_dir] - 1) )  ;
	    LatticeReal B = -Bfield*XX*Layout::lattSize()[y_dir]*
	      Layout::latticeCoordinate(y_dir)  ;
	    
	    U = cmplx(cos(B),sin(B));
	    return U ;
	  }
	}
	else//do nothing
	  QDPIO::cerr<<"OOOPS! don't know how to set up magnetic field in other than 4 dimenssions\n"; 
      }
      U = 1.0 ;
      QDPIO::cout<<__func__<<" A_"<<mu<<"= 0"<<endl ;
      END_CODE();
      return U;
    }


     // Construct linear Electric field. Needs fix up...
    LatticeComplex
    LinearElectricExternalField::operator()(int mu) const
    {
      START_CODE();
	  
      LatticeComplex U ;
      Real A0;
      A0 = -(0.5)*Efield;
      if(mu==t_dir){
	QDPIO::cout<<__func__<<" A_"<<mu<<"="<<A0<<endl ;
	// This line is -1/2 * Efield * (z-z0) * (z-z0-1) to get the correct E-field
	LatticeReal E = A0*(Layout::latticeCoordinate(x_dir)-x_src)*(Layout::latticeCoordinate(x_dir)-x_src - 1) ; 
	U = cmplx(cos(E),sin(E));
	return U ;
      }
      if((mu==x_dir)&&(patch)){// do the patch for the e-field
	LatticeBoolean XX = (Layout::latticeCoordinate(x_dir)==
			     (Layout::lattSize()[x_dir] - 1) )  ;
	LatticeReal E = -A0*XX*(Layout::lattSize()[x_dir]-1)*(Layout::lattSize()[x_dir]-2)*
	  Layout::latticeCoordinate(t_dir)  ;
	
	U = cmplx(cos(E),sin(E));
	return U ;
      }

      U = 1.0 ;
      QDPIO::cout<<__func__<<" A_"<<mu<<"= 0"<<endl ;
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
	//QDPIO::cerr<<"registered ZERO external field\n";

	success &= Chroma::TheExternalFieldFactory::Instance().registerObject(string("CONSTANT_MAGNETIC"), ConstantMagneticFunc);
	//QDPIO::cerr<<"registered CONSTANT_MAGNETIC field\n";

	success &= Chroma::TheExternalFieldFactory::Instance().registerObject(string("LINEAR_ELECTRIC"), LinearElectricFunc);
	//QDPIO::cerr<<"registered LINEAR_ELECTRIC field\n";

	registered = true;
      }
      return success;
    }
  }  // end namespace ExternalFieldEnv

}  // end namespace Chroma



  
