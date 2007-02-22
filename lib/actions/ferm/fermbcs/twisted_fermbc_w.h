// -*- C++ -*-
// $Id: twisted_fermbc_w.h,v 3.3 2007-02-22 21:11:45 bjoo Exp $
/*! \file
 *  \brief Twisted fermionic BC
 */

#ifndef __twisted_fermbc_w_h__
#define __twisted_fermbc_w_h__

#include "fermbc.h"
#include "handle.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma 
{
  //! Params struct for twisted params
  /*! \ingroup fermbcs */
  struct TwistedFermBCParams {
    TwistedFermBCParams() {}
    TwistedFermBCParams(XMLReader& in, const std::string& path);
    multi1d<int>  boundary;
    multi1d<Real> phases_by_pi;
    multi1d<int>  phases_dir;
  };
  
  // Readers writers
  void read(XMLReader& xml, const std::string& path, TwistedFermBCParams& param);
  void write(XMLWriter& xml, const std::string& path, const TwistedFermBCParams& param);


  //! Name and registration
  /*! \ingroup fermbcs */
  namespace WilsonTypeTwistedFermBCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Concrete class for all fermionic actions with twisted boundary conditions
  /*!
   * \ingroup fermbcs
   *
   *  Twisted BC
   */
  template<class T>
  class TwistedFermBC : public FermBC<T,
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

    TwistedFermBC(const multi1d<int>& boundary_,
		  const multi1d<Real>& phases_by_pi_,
		  const multi1d<int>& phases_dir_) : 
      phases_by_pi(phases_by_pi_),
      phases_dir(phases_dir_), 
      simple_bc_handle(new SimpleFermBC<T,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >(boundary_)) 
      {      
	check_arrays(phases_by_pi_, phases_dir_);
      }
   
    //! Copy constructor
    TwistedFermBC(const TwistedFermBC& a) :
      phases_by_pi(a.phases_by_pi),
      phases_dir(a.phases_dir), 
      simple_bc_handle(a.simple_bc_handle) 
      {
	check_arrays(a.phases_by_pi, a.phases_dir);    
      }
  
    //! Destructor is automatic
    ~TwistedFermBC() {}

    //! Assignment
    TwistedFermBC& operator=(const TwistedFermBC& a)
      { 
	phases_by_pi = a.phases_by_pi;
	phases_dir = a.phases_dir;
	simple_bc_handle = a.simple_bc_handle;
	check_arrays(a.phases_by_pi, a.phases_dir);    
	return *this;
      }


    //! Modify U fields in place
    void modify(multi1d<LatticeColorMatrix>& u) const
    {
      START_CODE();

      // Apply the simple BC's first
      (*simple_bc_handle).modify(u);
      const Real twopi = Real(6.283185307179586476925286);
      const Real onepi = twopi/Real(2);

      const multi1d<int>& nrow = Layout::lattSize();

      // Loop over the three twist angles
      for(int i=0; i < Nd-1; i++) 
      {
	int direction = phases_dir[i]; 
	Real arg = phases_by_pi[i]*onepi / Real( nrow[ direction ] );
	Real re = cos(arg);
	Real im = sin(arg);
	Complex phase = cmplx( re, im );
	u[ direction ] *= phase;
      }
    
      END_CODE();
    }

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(T& psi) const {}
 
    //! Modify fermion fields in place under a subset
    /*! NOP */
    void modifyF(T& psi, const Subset& s) const {}

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(multi1d<T>& psi) const {}
    
    //! Modify fermion fields in place under a subset
    /*! NOP */
    void modifyF(multi1d<T>& psi, const Subset& s) const {}

    //! Zero some gauge-like field in place on the masked links
    /*! NOP */
    void zero(multi1d<LatticeColorMatrix>& ds_u) const {}

    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return true;}


  private:
    // No empty constructor
    TwistedFermBC() {}

    void check_arrays( const multi1d<Real> phases_by_pi,
		       const multi1d<int> phases_dir ) 
    {
      START_CODE();
     
      if( phases_by_pi.size() != (Nd-1) ) { 
	QDPIO::cerr << "TwistedFermBCParams: Invalid size for phases_by_pi. Should be " << Nd-1 << "  but is " << phases_by_pi.size() << endl;
	QDP_abort(1);
      }
      
      if( phases_dir.size() != (Nd-1) ) { 
	QDPIO::cerr << "TwistedFermBCParams: Invalid size for phases_dir. Should be " << Nd-1 << "  but is " << phases_dir.size() << endl;
	QDP_abort(1);
      }
      
      for(int i=0; i < Nd-1; i++) { 
	if( toBool( phases_dir[i] < 0 || phases_dir[i] > Nd-1 ) ) { 
	  QDPIO::cerr << "Invalid value in phases_dir, direction " << i << " should be between 0 and " << Nd-1 << " but is " << phases_dir[i] << endl;
	}
      }
    
      END_CODE();
    }
   
  private:
    multi1d<Real> phases_by_pi;
    multi1d<int>  phases_dir;
    multi1d<int>  boundary;

    Handle< SimpleFermBC<T,
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> > > simple_bc_handle;
  };


} // Namespace Chroma 


// End of include guard 
#endif
