// -*- C++ -*-
// $Id: gaugebc_factory.h,v 3.1 2009-04-17 02:05:35 bjoo Exp $
/*! \file
 *  \brief Gauge boundary condition factories
 */

#ifndef __gaugebc_factory_h__
#define __gaugebc_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "gaugebc.h"

namespace Chroma
{

  //! GaugeAct Factory 
  /*! @ingroup gaugebcs */
  typedef Chroma::SingletonHolder< 
  ObjectFactory<GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    GaugeBC<multi1d<LatticeColorMatrix>, 
	    multi1d<LatticeColorMatrix> >* (*)(XMLReader&, const std::string&), 
		StringFactoryError> >
  TheGaugeBCFactory;

  typedef Chroma::SingletonHolder< 
  ObjectFactory<GaugeBC<multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    GaugeBC<multi1d<LatticeColorMatrixF>, 
	    multi1d<LatticeColorMatrixF> >* (*)(XMLReader&, const std::string&), 
		StringFactoryError> >
  TheGaugeBCFFactory;

  typedef Chroma::SingletonHolder< 
  ObjectFactory<GaugeBC<multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    GaugeBC<multi1d<LatticeColorMatrixD>, 
	    multi1d<LatticeColorMatrixD> >* (*)(XMLReader&, const std::string&), 
		StringFactoryError> >
  TheGaugeBCDFactory;
}; // end namespace Chroma


#endif
