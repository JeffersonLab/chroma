// -*- C++ -*-
/*! \file
 *  \brief Gamma insertions
 */

#ifndef __gamma_displacement_w_h__
#define __gamma_displacement_w_h__

#include "meas/smear/quark_displacement.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup smear */
  namespace GammaDisplacementEnv
  {
    bool registerAll();
  
    //! Return the name
    std::string getName();

    //! Params for simple quark displacement
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      int      gamma;          /*!< gamma value of insertion */
    };


    //! Gamma insertion/displacement
    /*! @ingroup smear
     *
     * Simple gamma multiplication of an object
     */
    template<typename T>
    class QuarkDisplace : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      QuarkDisplace(const Params& p) : params(p) {}
      
      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign) const
      {
	T tmp = Gamma(params.gamma) * quark;
	quark = tmp;
      }

    private:
      //! Hide partial constructor
      QuarkDisplace() {}

    private:
      Params  params;   /*!< displacement params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const std::string& path, GammaDisplacementEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const std::string& path, const GammaDisplacementEnv::Params& param);

}  // end namespace Chroma

#endif
