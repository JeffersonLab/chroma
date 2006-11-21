// -*- C++ -*-
// $Id: quark_flavor_s.h,v 1.2 2006-11-21 05:20:14 kostas Exp $
/*! \file
 *  \brief Staggered Flavor operators
 *
 *
 */

#ifndef __quark_flavor_s_h__
#define __quark_flavor_s_h__

#include "meas/smear/quark_displacement.h"
#include "util/ferm/staggered_operators_s.h"

namespace Chroma 
{
  using namespace StaggeredFlavorOperators;

  //! Name and registration
  /*! @ingroup smear */
  namespace StaggeredQuarkFlavorOpEnv
  {
    bool registerAll();

  
    //! Params for derivative quark displacement
    /*! @ingroup smear */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      FlavorOp;    /*! flavor operator type */
    };

    struct ParamsOneIndex
    {
      ParamsOneIndex();
      ParamsOneIndex(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      FlavorOp;    /*! flavor operator type */
      int     mu ; /*! flavor operator indices */
    };
    
    struct ParamsTwoIndex
    {
      ParamsTwoIndex();
      ParamsTwoIndex(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      FlavorOp;    /*! flavor operator type */
      int     mu ; /*! flavor operator indices */
      int     nu ; /*! flavor operator indices */
    };

    //! Construct staggered scalar flavored sources
    /*!
     * \ingroup sources
     *
     */
    template<typename T>
    class StaggeredScalarOp : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      StaggeredScalarOp(const Params& p) : params(p) {}
      
      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign ) const;

    private:
      Params  params;   /*!< source params */
    };

    //! Construct staggered pseudo scalar flavored sources
    /*!
     * \ingroup sources
     *
     */
    template<typename T>
    class StaggeredPseudoScalarOp : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      StaggeredPseudoScalarOp(const Params& p) : params(p) {}
      
      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign ) const;

    private:
      Params  params;   /*!< source params */
    };

    //! Construct staggered vector flavored sources
    /*!
     * \ingroup sources
     *
     */
    template<typename T>
    class StaggeredVectorOp : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      StaggeredVectorOp(const ParamsOneIndex& p) : params(p) {}
      
      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign ) const;

    private:
      ParamsOneIndex  params;   /*!< source params */
    };

    //! Construct staggered axial vector flavored sources
    /*!
     * \ingroup sources
     *
     */
    template<typename T>
    class StaggeredAxialVectorOp : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      StaggeredAxialVectorOp(const ParamsOneIndex& p) : params(p) {}
      
      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign ) const;

    private:
      ParamsOneIndex  params;   /*!< source params */
    };

    //! Construct tensor flavored sources
    /*!
     * \ingroup sources
     *
     */
    template<typename T>
    class StaggeredTensorOp : public QuarkDisplacement<T>
    {
    public:
      //! Full constructor
      StaggeredTensorOp(const ParamsTwoIndex& p) : params(p) {}
      
      //! Displace the quark
      void operator()(T& quark, 
		      const multi1d<LatticeColorMatrix>& u,
		      enum PlusMinus isign ) const;

    private:
      ParamsTwoIndex  params;   /*!< source params */
    };


  }  // end namespace


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, StaggeredQuarkFlavorOpEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const StaggeredQuarkFlavorOpEnv::Params& param);


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, StaggeredQuarkFlavorOpEnv::ParamsTwoIndex& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const StaggeredQuarkFlavorOpEnv::ParamsTwoIndex& param);


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, StaggeredQuarkFlavorOpEnv::ParamsOneIndex& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const StaggeredQuarkFlavorOpEnv::ParamsOneIndex& param);



}  // end namespace Chroma

#endif
