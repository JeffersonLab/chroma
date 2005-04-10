#ifndef SCHRGAUGEBC_PARAM_H
#define SCHRGAUGEBC_PARAM_H

#include "chromabase.h"
#include "gaugebc.h"
#include "io/enum_io/enum_gaugebc_io.h"



namespace Chroma { 
  

  /*! @ingroup gaugebcs */
  struct GaugeBCSchrParams {
    GaugeBCSchrParams();
    GaugeBCSchrParams(XMLReader& xml, const std::string& path);
    SchrFunType SchrFun;
    Real SchrPhiMult;
  };
  
  /*! @ingroup gaugebcs */
  void read(XMLReader& xml, const std::string& path, GaugeBCSchrParams& p);

  /*! @ingroup gaugebcs */
  void write(XMLWriter& xml, const std::string& path, const GaugeBCSchrParams& p) ;


  //! Concrete class for 1-link gauge action boundary conditions with Schroedinger BC
  /*! @ingroup gaugebcs
   *
   *  Schroedinger BC for gauge actions that have only 1 link for padding
   *  in the decay direction - e.g., Wilson-like gauge action
   */
  class Schr1LinkGaugeBC : public SchrGaugeBC
  {
  public:
    //! Schroedinger BC
    /*! 
     * \param SchrFun_      type of Schroedinger BC
     * \param SchrPhiMult_  factor to rescale fixed field
     */
    Schr1LinkGaugeBC(const SchrFunType SchrFun_, const Real& SchrPhiMult_)
      {
	initFunc(SchrFun_, SchrPhiMult_);
      }



    //! Initialise from params struct
    Schr1LinkGaugeBC(const GaugeBCSchrParams& params) 
      {
	initFunc(params.SchrFun, params.SchrPhiMult);
      }


    //! Copy constructor
    Schr1LinkGaugeBC(const Schr1LinkGaugeBC& a) : SchrFun(a.SchrFun), 
      decay_dir(a.decay_dir), mask(a.mask), fld(a.fld) {}

    //! Destructor is automatic
    ~Schr1LinkGaugeBC() {}

    //! Assignment
    Schr1LinkGaugeBC& operator=(const Schr1LinkGaugeBC& a)
      {
	SchrFun = a.SchrFun; decay_dir = a.decay_dir; 
	mask = a.mask; fld = a.fld;
	return *this;
      }

#if defined(EXPOSE_THIS_STUFF)
    //! Type of Schroedinger BC
    SchrFunType getSFBC() const {return SchrFun;}
#endif

    //! Modify U fields in place
    void modify(multi1d<LatticeColorMatrix>& u) const
      {QDP_error_exit("modify not implemented");}

    //! Zero the U fields in place on the masked links
    void zero(multi1d<LatticeColorMatrix>& u) const
      {QDP_error_exit("zero not implemented");}

#if defined(EXPOSE_THIS_STUFF)
    // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lbmaskU() const {return mask;}

    //! Fixed gauge links on only the lbmaskU() sites
    const multi1d<LatticeColorMatrix>& lFldU() const {return fld;}
#endif

    //! Says if there are fixed links within the lattice
    bool nontrivialP() const {return true;}

    //! Decay direction
    int getDir() const {return decay_dir;}

  private:
    // Hide default constuctor
    Schr1LinkGaugeBC() {}
    void initFunc(const SchrFunType SchrFun_, const Real& SchrPhiMult_) {
      QDPIO::cerr << "Schr1LinkGaugeBC() not yet implemented" << endl;
      QDP_abort(1);
    }

  private:
    SchrFunType SchrFun;
    int decay_dir;
    multi1d<LatticeBoolean> mask;
    multi1d<LatticeColorMatrix> fld;
  };


  //! Concrete class for 2-link gauge action boundary conditions with Schroedinger BC
  /*! @ingroup actions
   *
   *  Schroedinger BC for gauge actions that have only 2 links for padding
   *  in the decay direction - e.g., Symanzik-like gauge action
   */
  class Schr2LinkGaugeBC : public SchrGaugeBC
  {
  public:
    //! Schroedinger BC
    /*! 
     * \param SchrFun_      type of Schroedinger BC
     * \param SchrPhiMult_  factor to rescale fixed field
     */
    Schr2LinkGaugeBC(const SchrFunType SchrFun_, const Real& SchrPhiMult_)
      {
	initFunc(SchrFun_, SchrPhiMult_);
      }

    //! Initialise from params struct
    Schr2LinkGaugeBC(const GaugeBCSchrParams& params)
      {
	initFunc(params.SchrFun, params.SchrPhiMult);
      }

    //! Copy constructor
    Schr2LinkGaugeBC(const Schr2LinkGaugeBC& a) : SchrFun(a.SchrFun), 
      decay_dir(a.decay_dir), mask(a.mask), fld(a.fld) {}

    //! Destructor is automatic
    ~Schr2LinkGaugeBC() {}

    //! Assignment
    Schr2LinkGaugeBC& operator=(const Schr2LinkGaugeBC& a)
      {
	SchrFun = a.SchrFun; decay_dir = a.decay_dir; 
	mask = a.mask; fld = a.fld;
	return *this;
      }

#if defined(EXPOSE_THIS_STUFF)
    //! Type of Schroedinger BC
    SchrFunType getSFBC() const {return SchrFun;}
#endif

    //! Modify U fields in place
    void modify(multi1d<LatticeColorMatrix>& u) const
      {QDP_error_exit("modify not implemented");}

    //! Zero the U fields in place on the masked links
    void zero(multi1d<LatticeColorMatrix>& u) const
      {QDP_error_exit("zero not implemented");}

#if defined(EXPOSE_THIS_STUFF)
    // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lbmaskU() const {return mask;}

    //! Fixed gauge links on only the lbmaskU() sites
    const multi1d<LatticeColorMatrix>& lFldU() const {return fld;}
#endif

    //! Says if there are fixed links within the lattice
    bool nontrivialP() const {return true;}

    //! Decay direction
    int getDir() const {return decay_dir;}

  private:
    // Hide default constuctor
    Schr2LinkGaugeBC() {}

    void initFunc(const SchrFunType SchrFun_, const Real& SchrPhiMult_) {
      QDPIO::cerr << "Schr2LinkGaugeBC() not yet implemented" << endl;
      QDP_abort(1);
    }

  private:
    SchrFunType SchrFun;
    int decay_dir;
    multi1d<LatticeBoolean> mask;
    multi1d<LatticeColorMatrix> fld;
  };

  
};

#endif
