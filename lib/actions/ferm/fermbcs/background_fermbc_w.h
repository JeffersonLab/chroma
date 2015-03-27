// -*- C++ -*-
/*! \file
 *  \brief quark bilinear background field fermionic BC
 */

#ifndef __background_fermbc_w_h__
#define __background_fermbc_w_h__

#include "fermbc.h"
#include "handle.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma 
{
  //! Params struct for twisted params
  /*! \ingroup fermbcs */
  struct BackgroundFermBCParams {
    BackgroundFermBCParams() {}
    BackgroundFermBCParams(XMLReader& in, const std::string& path);
    multi1d<int>  boundary;
    int gamma ; //the gamma matrix in the quark bi-linear
    Complex lambda; // the field
  };
  
  // Readers writers
  void read(XMLReader& xml, const std::string& path, BackgroundFermBCParams& param);
  void write(XMLWriter& xml, const std::string& path, const BackgroundFermBCParams& param);


  //! Name and registration
  /*! \ingroup fermbcs */
  namespace WilsonTypeBackgroundFermBCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Concrete class for all fermionic actions with twisted boundary conditions
  /*!
   * \ingroup fermbcs
   *
   *  Background BC
   */
  template<class T>
  class BackgroundFermBC : public FermBC<T,
				      multi1d<LatticeColorMatrix>, 
				      multi1d<LatticeColorMatrix> >
  {
  public:

    //! Only full constructor
    /*!
     * \param boundary_phases  multiply links on edge of lattice by boundary
     *
     * NOTE: there is no real reason this is of type int, could be more general
     *       like Complex
     */

    BackgroundFermBC(int gamma_,
		     Complex lambda_,
		     const multi1d<int>& boundary_) : 
      gamma(gamma_),
      lambda(lambda_), 
      simple_bc_handle(new SimpleFermBC<T,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >(boundary_)) 
    {           }
    
    //! Copy constructor
    BackgroundFermBC(const BackgroundFermBC& a) :
      gamma(a.gamma),
      lambda(a.lambda), 
      simple_bc_handle(a.simple_bc_handle) 
      {      }
  
    //! Destructor is automatic
    ~BackgroundFermBC() {}

    //! Assignment
    BackgroundFermBC& operator=(const BackgroundFermBC& a)
      { 
	lambda = a.lambda;
	gamma = a.gamma ;
	simple_bc_handle = a.simple_bc_handle;
	return *this;
      }


    //! Modify U fields in place
    void modify(multi1d<LatticeColorMatrix>& u) const
    {
      START_CODE();

      // Apply the simple BC's 
      (*simple_bc_handle).modify(u);
        
      //nothing else to do here

      END_CODE();
    }

    //! Modify fermion fields in place
    /*! add the quark bilinear */
    void modifyF(T& psi) const {
      psi +=  lambda*(Gamma(gamma)*psi) ;
    }
 
    //! Modify fermion fields in place under a subset
    /*! add the quark bilinear  */
    void modifyF(T& psi, const Subset& s) const {
      psi[s] +=  lambda*(Gamma(gamma)*psi) ;
    }

    //! Modify fermion fields in place
    /*! add the quark bilinear */
    void modifyF(multi1d<T>& psi) const {
      QDPIO::cout<<"BACKGROUND BC NOT IMPLEMENTED"<<std::endl ;
      QDP_abort(12003);
    }
    
    //! Modify fermion fields in place under a subset
    /*! add the quark bilinear */
    void modifyF(multi1d<T>& psi, const Subset& s) const {
      QDPIO::cout<<"BACKGROUND BC NOT IMPLEMENTED"<<std::endl ;
      QDP_abort(12003);
    }

    //! Zero some gauge-like field in place on the masked links
    /*! NOP */
    void zero(multi1d<LatticeColorMatrix>& ds_u) const {}

    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return true;}


  private:
    // No empty constructor
    BackgroundFermBC() {}

   
  private:
    int gamma;
    Complex lambda ;

    Handle< SimpleFermBC<T,
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> > > simple_bc_handle;
  };


} // Namespace Chroma 


// End of include guard 
#endif
