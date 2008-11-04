// -*- C++ -*-
// $Id: simple_quark_displacement.h,v 3.2 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Simple quark displacement
 */

#ifndef __simple_quark_displacement_h__
#define __simple_quark_displacement_h__

#include "meas/smear/quark_displacement.h"
#include "meas/smear/displacement.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup smear */
  namespace SimpleQuarkDisplacementEnv
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
    
      int              disp_length;          /*!< displacement length */
      int              disp_dir;             /*!< x(0), y(1), z(2) */   
    };


    //! Simple quark displacement
    /*! @ingroup smear
     *
     * Simple quark displacement object
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
	switch (isign)
	{
	case PLUS:
	  displacement(u, quark, params.disp_length, params.disp_dir);
	  break;

	case MINUS:
	  displacement(u, quark, -params.disp_length, params.disp_dir);
	  break;

	default:
	  QDP_error_exit("illegal isign in QuarkDisplace");
	}
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
  void read(XMLReader& xml, const string& path, SimpleQuarkDisplacementEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const SimpleQuarkDisplacementEnv::Params& param);

}  // end namespace Chroma

#endif
