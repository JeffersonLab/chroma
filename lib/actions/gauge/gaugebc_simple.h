#ifndef GAUGEBC_SIMPLE_H
#define GAUGEBC_SIMPLE_H

#include "chromabase.h"
#include "gaugebc.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

namespace Chroma { 
  
  namespace SimpleGaugeBCEnv { 
    extern const std::string name;
    extern const bool registered;
  }

  struct SimpleGaugeBCParams { 
    SimpleGaugeBCParams();
    SimpleGaugeBCParams(XMLReader& xml, const std::string& path);
    multi1d<Complex> boundary;
  };

  void read(XMLReader& xml, const std::string& path, SimpleGaugeBCParams& p); 
  
  void write(XMLWriter& xml, const std::string& path, const SimpleGaugeBCParams& p);

 //! Concrete class for gauge actions with simple boundary conditions
  /*! @ingroup actions
   *
   *  Simple BC, where boundary array is multiplied on the links on the
   *  edge of the lattice
   */
  class SimpleGaugeBC : public GaugeBC
  {
  public:
    //! Only full constructor
    /*!
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * NOTE: the type of boundary can be generalized
     */
    SimpleGaugeBC(const multi1d<Complex>& boundary_) : boundary(boundary_) {}

    //! From param struct 
    SimpleGaugeBC(const SimpleGaugeBCParams& p) : boundary(p.boundary) {}

    //! Copy constructor
    SimpleGaugeBC(const SimpleGaugeBC& a) : boundary(a.boundary) {}

    //! Destructor is automatic
    ~SimpleGaugeBC() {}

    //! Assignment
    SimpleGaugeBC& operator=(const SimpleGaugeBC& a)
      {boundary = a.boundary; return *this;}

    //! Modify U fields in place
    void modify(multi1d<LatticeColorMatrix>& u) const
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
    void zero(multi1d<LatticeColorMatrix>& u) const {}

    //! Says if there are non-trivial BC links
    /*! 
     * This is a simple implementation - always do the work. 
     * Could be improved by checking boundary
     */
    bool nontrivialP() const {return true;}

  private:
    // Hide empty constructor
    SimpleGaugeBC() {}

  private:
    multi1d<Complex> boundary;
  };

};


#endif
