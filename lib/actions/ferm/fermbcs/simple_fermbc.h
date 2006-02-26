// -*- C++ -*-
// $Id: simple_fermbc.h,v 2.2 2006-02-26 03:47:52 edwards Exp $
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
  /*! \ingroup fermbc */
  struct SimpleFermBCParams
  {
    SimpleFermBCParams() {}
    SimpleFermBCParams(XMLReader& in, const std::string& path);

    multi1d<int> boundary;
  };

  // Reader/writers
  /*! \ingroup fermbc */
  void read(XMLReader& xml, const std::string& path, SimpleFermBCParams& param);
  /*! \ingroup fermbc */
  void write(XMLWriter& xml, const std::string& path, const SimpleFermBCParams& param);


  //! Concrete class for all gauge actions with simple boundary conditions
  /*! @ingroup fermbc
   *
   *  Simple BC, where boundary array is multiplied on the links on the
   *  edge of the lattice
   */
  template<class T>
  class SimpleFermBC : public FermBC<T>
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
      gbc(new PeriodicGaugeBC()), boundary(p.boundary)
      {init();}
  
    //! Only full constructor
    /*!
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * NOTE: there is no real reason this is of type int, could be more general
     *       like Complex
     */
    SimpleFermBC(const multi1d<int>& boundary_) : 
      gbc(new PeriodicGaugeBC()), boundary(boundary_)
      {init();}
  
    //! Use both a Gauge BC and a boundary field
    /*!
     * \param gbc       A simple Gauge BC
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * The type of gbc is restricted here mainly for paranoia.
     * The type could be loosened (in fact to an abstract type), but it makes
     * no sense to apply these type FermBC to a Schr1LinkGaugeBC, for example.
     *
     * Currently, no code uses this capability...
     */
    SimpleFermBC(const SimpleGaugeBC& gbc_, const multi1d<int>& boundary_) : 
      gbc(new SimpleGaugeBC(gbc_)), boundary(boundary_)
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
    void modifyU(multi1d<LatticeColorMatrix>& u) const
      {
	gbc->modify(u);   // first modify the U fields

	if (nontrivialP())
	{
	  // phases for all the directions
	  for(int m = 0; m < Nd; ++m)
	  {
	    /* u[m] *=  (coord(m) == nrow(m)-1 ) ? boundary[m] : 1 */
	    u[m] *= where(Layout::latticeCoordinate(m) == Layout::lattSize()[m]-1,
			  Integer(boundary[m]), Integer(1));
	  }
	}
      }

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(T& psi) const {}
 
    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return nontriv;}

  protected:
    void init()
      {
	// Check triviality of boundary
	nontriv = false;
	for(int m = 0; m < Nd; ++m)
	  if (boundary[m] != 1)
	    nontriv |= true;
      }

  private:
    // No empty constructor
    SimpleFermBC() {}

  private:
    Handle<GaugeBC> gbc;
    multi1d<int> boundary;
    bool nontriv;
  };

}

#endif
