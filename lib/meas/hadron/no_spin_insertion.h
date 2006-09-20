// -*- C++ -*-
// $Id: no_spin_insertion.h,v 1.2 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief No spin insertion
 */

#ifndef __no_spin_insertion_h__
#define __no_spin_insertion_h__

#include "meas/hadron/spin_insertion.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup hadron */
  namespace NoSpinInsertionEnv
  {
    extern const std::string name;
    bool registerAll();
  

    //! Params for no spin insertion
    /*! @ingroup hadron */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    };


    //! No spin insertion
    /*! @ingroup hadron
     *
     * No spin insertion object
     */
    template<typename T>
    class SpinInsert : public SpinInsertion<T>
    {
    public:
      //! Full constructor
      SpinInsert(const Params& p) : params(p) {}
      
      //! Displace the quark
      T operator()(const T& quark) const {return quark;}

    private:
      //! Hide partial constructor
      SpinInsert() {}

    private:
      Params  params;   /*!< displacement params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, NoSpinInsertionEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const NoSpinInsertionEnv::Params& param);

}  // end namespace Chroma

#endif
