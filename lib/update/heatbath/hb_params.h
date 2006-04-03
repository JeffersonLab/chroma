// -*- C++ -*-
// $Id: hb_params.h,v 3.0 2006-04-03 04:59:07 edwards Exp $

/*! \file
 * \brief Parameters for heat-bath
 */

#ifndef HB_PARAMS_H
#define HB_PARAMS_H

namespace Chroma 
{

  //! Heat-bath params
  /*! \ingroup heatbath */
  struct HBParams 
  {
    int nmax() const { return NmaxHB; }
    Double beta() const { return BetaMC; }
    Double xi() const { return xi_0; }
    Double xi2() const { return xi_0*xi_0; }
    bool aniso() const {return anisoP; }

    /**************************************************
     * number of maximum HB tries for Creutz or KP a_0, 
     * negative or zero value - update every single link
     * (try infinitely long)
     **************************************************/
    int NmaxHB;
    // MC SU(N) Beta
    Double BetaMC;
    //the bare anisotropy
    Double xi_0;
    int  t_dir;
    int  nOver;
    bool anisoP;
  };

  
  //! Reader
  void read(XMLReader& xml, const std::string& path, HBParams& p);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const HBParams& p);

}  // end namespace Chroma

#endif
