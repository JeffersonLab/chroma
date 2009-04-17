// -*- C++ -*-
// $Id: simple_fermbc.h,v 3.4 2009-04-17 02:05:30 bjoo Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_h__
#define __simple_fermbc_h__

#include "fermbc.h"
#include "handle.h"
#include "actions/gauge/gaugebcs/simple_gaugebc.h"
#include "actions/gauge/gaugebcs/periodic_gaugebc.h"

namespace Chroma
{

  //! Params for simple fermbc
  /*! \ingroup fermbcs */
  struct SimpleFermBCParams
  {
    SimpleFermBCParams() {}
    SimpleFermBCParams(XMLReader& in, const std::string& path);

    multi1d<int> boundary;
  };

  // Reader/writers
  /*! \ingroup fermbcs */
  void read(XMLReader& xml, const std::string& path, SimpleFermBCParams& param);
  /*! \ingroup fermbcs */
  void write(XMLWriter& xml, const std::string& path, const SimpleFermBCParams& param);


  //! Concrete class for all gauge actions with simple boundary conditions
  /*! @ingroup fermbcs
   *
   *  Simple BC, where boundary array is multiplied on the links on the
   *  edge of the lattice
   */
  template<class T, typename P, typename Q>
  class SimpleFermBC : public FermBC<T,P,Q>
  {
  public:
    //! Only full constructor
    /*!
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * NOTE: there is no real reason this is of type int, could be more general
     *       like Complex
     */
    SimpleFermBC(const SimpleFermBCParams& p) :
      gbc(new PeriodicGaugeBC<P,Q>()), boundary(p.boundary)
      {init();}
  
    //! Only full constructor
    /*!
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * NOTE: there is no real reason this is of type int, could be more general
     *       like Complex
     */
    SimpleFermBC(const multi1d<int>& boundary_) : 
      gbc(new PeriodicGaugeBC<P,Q>()), boundary(boundary_)
      {init();}
  
    //! Use both a Gauge BC and a boundary field
    /*!
     * \param gbc       A generic gaugebc
     * \param boundary  multiply links on edge of lattice by boundary
     */
    SimpleFermBC(const Handle< GaugeBC<P,Q> >& gbc_, const multi1d<int>& boundary_) : 
      gbc(gbc_), boundary(boundary_)
      {init();}
  
    //! Copy constructor
    SimpleFermBC(const SimpleFermBC& a) : 
      gbc(a.gbc), boundary(a.boundary), nontriv(a.nontriv) {}

    //! Destructor is automatic
    ~SimpleFermBC() {}

    //! Assignment
    SimpleFermBC& operator=(const SimpleFermBC& a)
      {gbc = a.gbc; boundary = a.boundary; nontriv = a.nontriv; return *this;}

    //! Modify U fields in place
    void modify(Q& u) const
    {
      START_CODE();

      gbc->modify(u);   // first modify the U fields

      // phases for all the directions
      for(int m = 0; m < Nd; ++m) 
      { 
	if( boundary[m] != 1 ) { 
	  /* u[m] *=  (coord(m) == nrow(m)-1 ) ? boundary[m] : 1 */
	  u[m] *= where(Layout::latticeCoordinate(m) == (Layout::lattSize()[m]-1),
			Integer(boundary[m]), Integer(1));
	}
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
    void zero(P& ds_u) const {}

    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return false;}

  protected:
    void init() {}

  private:
    // No empty constructor
    SimpleFermBC() {}

  private:
    Handle< GaugeBC<P,Q> > gbc;
    multi1d<int> boundary;
    bool nontriv;
  };

}

#endif
