// -*- C++ -*-
// $Id: simple_spin_insertion_w.h,v 1.2 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief Gamma insertions
 */

#ifndef __simple_gamma_insertion_w_h__
#define __simple_gamma_insertion_w_h__

#include "meas/hadron/spin_insertion.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleSpinInsertionEnv
  {
    extern const std::string name;
    bool registerAll();
  

    //! Params for simple spin insertion
    /*! @ingroup hadron */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      int      gamma;          /*!< gamma value of insertion */
    };


    //! Gamma insertion
    /*! @ingroup hadron
     *
     * Simple gamma multiplication of an object
     */
    template<typename T>
    class LeftSpinInsert : public SpinInsertion<T>
    {
    public:
      //! Full constructor
      LeftSpinInsert(const Params& p) : params(p) {}
      
      //! Spin insert
      T operator()(const T& quark) const {return Gamma(params.gamma) * quark;}

    private:
      //! Hide partial constructor
      LeftSpinInsert() {}

    private:
      Params  params;   /*!< spin insertion params */
    };


    //! Gamma insertion
    /*! @ingroup hadron
     *
     * Simple gamma multiplication of an object
     */
    template<typename T>
    class RightSpinInsert : public SpinInsertion<T>
    {
    public:
      //! Full constructor
      RightSpinInsert(const Params& p) : params(p) {}
      
      //! Spin insert
      T operator()(const T& quark) const {return quark * Gamma(params.gamma);}

    private:
      //! Hide partial constructor
      RightSpinInsert() {}

    private:
      Params  params;   /*!< spin insertion params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleSpinInsertionEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const SimpleSpinInsertionEnv::Params& param);

}  // end namespace Chroma

#endif
