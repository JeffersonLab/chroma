// -*- C++ -*-
/*! \file
 *  \brief No quark displacement
 */

#ifndef __no_quark_displacement_h__
#define __no_quark_displacement_h__

#include "meas/smear/quark_displacement.h"
#include "meas/smear/displacement.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup smear */
  namespace NoQuarkDisplacementEnv
  {
    bool registerAll();
  
    //! Return the name
    std::string getName();

    //! Params for no quark displacement
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    };


    //! No quark displacement
    /*! @ingroup smear
     *
     * No quark displacement object
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
	{ /* nop */ }

    private:
      //! Hide partial constructor
      QuarkDisplace() {}

    private:
      Params  params;   /*!< displacement params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const std::string& path, NoQuarkDisplacementEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const std::string& path, const NoQuarkDisplacementEnv::Params& param);

}  // end namespace Chroma

#endif
