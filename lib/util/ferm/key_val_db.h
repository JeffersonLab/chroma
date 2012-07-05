// -*- C++ -*-
/*! \file
 * \brief Key and values for DB
 */

#ifndef __key_val_db_h__
#define __key_val_db_h__

#include "chromabase.h"
#include "qdp_db.h"

namespace Chroma
{
  using namespace FILEDB;

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

    void writeObject (std::string& output) const throw (SerializeException) {
      BinaryBufferWriter bin;
      write(bin, key());
      output = bin.strPrimaryNode();
    }

    void readObject (const std::string& input) throw (SerializeException) {
      BinaryBufferReader bin(input);
      read(bin, key());
    }

    // Part of DBKey
    int hasHashFunc (void) const {return 0;}
    int hasCompareFunc (void) const {return 0;}

    /**
     * Empty hash and compare functions. We are using default functions.
     */
    static unsigned int hash (const void* bytes, unsigned int len) {return 0;}
    static int compare (const FFDB_DBT* k1, const FFDB_DBT* k2) {return 0;}
   
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

    void writeObject (std::string& output) const throw (SerializeException) {
      BinaryBufferWriter bin;
      write(bin, data());
      output = bin.strPrimaryNode();
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
