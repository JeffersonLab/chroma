// -*- C++ -*-
// $Id: key_val_db.h,v 1.1 2008-08-05 04:15:04 edwards Exp $
/*! \file
 * \brief Key and values for DB
 */

#ifndef __key_val_db_h__
#define __key_val_db_h__

#include "chromabase.h"
#include "qdp_db.h"

namespace Chroma
{
  using namespace FFDB;

  //---------------------------------------------------------------------
  //! Serializable key harness
  /*! \ingroup ferm */
  template<typename K>
  class SerialDBKey : public DBKey
  {
  public:
    //! Default constructor
    SerialDBKey() {} 

    //! Constructor from data
    SerialDBKey(const K& k) : key_(k) {}

    //! Setter
    K& key() {return key_;}

    //! Getter
    const K& key() const {return key_;}

    // Part of Serializable
    const unsigned short serialID (void) const {return 456;}

    void writeObject (std::string& output) throw (SerializeException) {
      BinaryBufferWriter bin;
      write(bin, key());
      output = bin.str();
    }

    void readObject (const std::string& input) throw (SerializeException) {
      BinaryBufferReader bin(input);
      read(bin, key());
    }

    // Part of DBKey
    int hasHashFunc (void) const {return 0;}
    int hasCompareFunc (void) const {return 0;}

    /**
     * Static Hash Function Implementation
     */
    static unsigned int hash (Db *db, const void* bytes, unsigned int len) {return 0;}

    /**
     * Static empty compare function 
     */
    static int compare (Db *db, const Dbt* k1, const Dbt* k2) {return 0;}

  private:
    K  key_;
  };


  //---------------------------------------------------------------------
  //! Serializable value harness
  /*! \ingroup ferm */
  template<typename D>
  class SerialDBData : public DBData
  {
  public:
    //! Default constructor
    SerialDBData() {} 

    //! Constructor from data
    SerialDBData(const D& d) : data_(d) {}

    //! Setter
    D& data() {return data_;}

    //! Getter
    const D& data() const {return data_;}

    // Part of Serializable
    const unsigned short serialID (void) const {return 123;}

    void writeObject (std::string& output) throw (SerializeException) {
      BinaryBufferWriter bin;
      write(bin, data());
      output = bin.str();
    }

    void readObject (const std::string& input) throw (SerializeException) {
      BinaryBufferReader bin(input);
      read(bin, data());
    }

  private:
    D  data_;
  };

} // namespace Chroma

#endif
