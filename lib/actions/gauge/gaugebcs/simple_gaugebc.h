// -*- C++ -*-
// $Id: simple_gaugebc.h,v 3.3 2009-04-17 02:05:35 bjoo Exp $
/*! \file
 *  \brief Simple gauge boundary conditions
 */

#ifndef __simple_gaugebc_h__
#define __simple_gaugebc_h__

#include "gaugebc.h"

namespace Chroma 
{ 
  
  /*! @ingroup gaugebcs */
  namespace SimpleGaugeBCEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

  /*! @ingroup gaugebcs */
  struct SimpleGaugeBCParams 
  { 
    SimpleGaugeBCParams();
    SimpleGaugeBCParams(XMLReader& xml, const std::string& path);
    multi1d<Complex> boundary;
  };

  /*! @ingroup gaugebcs */
  void read(XMLReader& xml, const std::string& path, SimpleGaugeBCParams& p); 
  
  /*! @ingroup gaugebcs */
  void write(XMLWriter& xml, const std::string& path, const SimpleGaugeBCParams& p);


  //! Concrete class for gauge actions with simple boundary conditions
  /*! @ingroup gaugebcs
   *
   *  Simple BC, where boundary array is multiplied on the links on the
   *  edge of the lattice
   */
  template<typename P, typename Q>
  class SimpleGaugeBC : public GaugeBC<P,Q>
  {
  public:
    // Typedefs to save typing
    //! Only full constructor
    /*!
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * NOTE: the type of boundary can be generalized
     */
    SimpleGaugeBC(const multi1d<Complex>& boundary_) : boundary(boundary_) {}

    //! From param struct 
    SimpleGaugeBC(const SimpleGaugeBCParams& p) : boundary(p.boundary) {}

    //! Destructor is automatic
    ~SimpleGaugeBC() {}

    //! Modify U fields in place
    void modify(Q& u) const
      {
	// phases for all the directions
	for(int m = 0; m < Nd; ++m)
	{
	  /* u[m] *=  (coord(m) == nrow(m)-1 ) ? boundary[m] : 1 */
	  u[m] *= where(Layout::latticeCoordinate(m) == Layout::lattSize()[m]-1,
			boundary[m], Complex(1));
	}
      }

    //! Zero the U fields in place on the masked links
    void zero(P& u) const {}

    //! Says if there are non-trivial BC links
    /*! 
     * This is a simple implementation - always do the work. 
     * Could be improved by checking boundary
     */
    bool nontrivialP() const {return false;}

  private:
    // Hide empty constructor
    SimpleGaugeBC() {}

    //! Hide assignment
    void operator=(const SimpleGaugeBC<P,Q>& a) {}

  private:
    multi1d<Complex> boundary;
  };

};


#endif
