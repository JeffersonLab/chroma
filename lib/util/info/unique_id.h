// -*- C++ -*-
/*! \file
 *  \brief Generate a unique id
 */

#ifndef __unique_id_h__
#define __unique_id_h__

namespace Chroma 
{

  //! Return a unique id
  /*!
   * \ingroup info
   *
   *  The id is return type a std::string. This gives maximal flexibility allowing
   *  the way the ID is generated to change in the future.
   */
  std::string uniqueId();

}  // end namespace Chroma

#endif
