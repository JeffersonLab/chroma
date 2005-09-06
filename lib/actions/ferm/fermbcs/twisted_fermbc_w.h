// Include guard
#ifndef twisted_fermbc_w_h
#define twisted_fermbc_w_h

#include "fermbcs.h"

namespace Chroma { 

  //! Params struct for twisted params
  struct TwistedFermBCParams {
    TwistedFermBCParams() {}
    TwistedFermBCParams(XMLReader& in, const std::string& path);
    multi1d<Real> boundary_phases;
  };

  namespace TwistedFermBCEnv {
    extern const std::string name;
  }

  // Readers writers
  void read(XMLReader& xml, const std::string& path, TwistedFermBCParams& param);
  void write(XMLWriter& xml, const std::string& path, const TwistedFermBCParams& param);


  // 4d & 5D part
  //! 4d name and registration
  namespace WilsonTypeTwistedFermBCEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! 5d Name and registration
  namespace WilsonTypeTwistedFermBCArrayEnv {
    extern const std::string name;
    extern const bool registered;
  }

  //! Concrete class for all gauge actions with twisted boundary conditions
  /*! @ingroup actions
   *
   *  Twisted BC, 
   */
  template<class T>
  class TwistedFermBC : public FermBC<T>
  {
  public:
    //! Only full constructor
    /*!
     * \param boundary_phases  multiply links on edge of lattice by boundary
     *
     * NOTE: there is no real reason this is of type int, could be more general
     *       like Complex
     */

    TwistedFermBC(const multi1d<Real>& boundary_phases_) : 
      boundary_phases(boundary_phases_) {
      QDPIO::cout << "TWISTED BOUNDARIES BEING CREATED" << endl;

    }
   
    //! Copy constructor
    TwistedFermBC(const TwistedFermBC& a) : 
      boundary_phases(a.boundary_phases) {}

    //! Destructor is automatic
    ~TwistedFermBC() {}

    //! Assignment
    TwistedFermBC& operator=(const TwistedFermBC& a)
    { 
      boundary_phases = a.boundary_phases; 
      return *this;
    }


    //! Modify U fields in place
    void modifyU(multi1d<LatticeColorMatrix>& u) const
      {
	// Modify the U field here.
      }

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(T& psi) const {}
 
    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return true;}


  private:
    // No empty constructor
    TwistedFermBC() {}

  private:
    multi1d<Real> boundary_phases;
  };

}; // Namespace Chroma 


// End of include guard 
#endif
