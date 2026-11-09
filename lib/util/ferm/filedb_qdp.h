// -*- C++ -*-
/*! @file
 * @brief Support for filedb
 */

#ifndef CHROMA_FILEDB_QDP_H
#define CHROMA_FILEDB_QDP_H

#include "chromabase.h"
#include "ConfDataStoreDB.h"

namespace Chroma
{
  /*! @defgroup io IO
   *
   * FileDB support
   *
   * @{
   */

  /**
   * Database cache size in bytes
   */
  const int db_cachesize = 50*1024*1024;

  /**
   * Database page size in bytes
   */
  const int db_pagesize = 64*1024;

  //--------------------------------------------------------------------------------
  //!  Local DB Base class
  /*!
    This class is used for writing of user data (most usefully measurements)
    into a DB file with a key/value semantics.
    NOTE: NOT COLLECTIVE, all operations done on the scope of the current process.
          Use BinaryStoreDB for a collective behaviour.
  */
  template<typename K, typename D>
  class LocalBinaryStoreDB
  {
  public:
    /**
     * Constructor
     *
     * @param do_init: whether to set cache
     */
    LocalBinaryStoreDB(bool do_init = false)
    {
      if (do_init)
      {
	db.setCacheSize(db_cachesize);
	// db.enablePageMove();   // Disable!
	db.setMaxNumberConfigs(1);
      }
    }

    /*!
      Destroy the object
    */
    ~LocalBinaryStoreDB()
    {
      close();
    }

    /**
     * How much data and keys should be kept in memory in bytes
     *
     * This should be called before the open is called
     * @param max_cache_size number of bytes of data and keys should be kept
     * in memory
     */
    virtual void setCacheSize (const unsigned int size)
    {
      db.setCacheSize(size);
    }


    /**
     * How much data and keys should be kept in memory in megabytes
     *
     * This should be called before the open is called
     * @param max_cache_size number of bytes of data and keys should be kept
     * in memory
     */
    virtual void setCacheSizeMB (const unsigned int size)
    {
      db.setCacheSizeMB(size);
    }


    /**
     * Page size used when a new data based is created
     * This only effects a new database
     *
     * @param pagesize the pagesize used in hash database. This value
     * should be power of 2. The minimum value is 512 and maximum value
     * is 262144
     *
     */
    virtual void setPageSize (const unsigned int size)
    {
      db.setPageSize(size);
    }


    /**
     * Set initial number of buckets
     *
     * This should be called before the open is called
     *
     * @param num the number of buckets should be used
     *
     */
    virtual void setNumberBuckets (const unsigned int num)
    {
      db.setNumberBuckets(num);
    }


    /**
     * Set whether to move pages when close to save disk space
     *
     * This only effective on writable database
     *
     */
    virtual void enablePageMove (void)
    {
      db.enablePageMove();
    }

    virtual void disablePageMove (void)
    {
      db.disablePageMove();
    }

    /**
     * Set and get maximum user information length
     */
    virtual void setMaxUserInfoLen (unsigned int len)
    {
      db.setMaxUserInfoLen(len);
    }

    virtual unsigned int getMaxUserInfoLen (void) const
    {
      return db.getMaxUserInfoLen();
    }

    /**
     * Set and get maximum number of configurations
     */
    virtual void setMaxNumberConfigs (unsigned int num)
    {
      db.setMaxNumberConfigs(num);
    }

    virtual unsigned int getMaxNumberConfigs (void) const
    {
      return db.getMaxNumberConfigs();
    }


    /**
     * Check if a DB file exists before opening.
     */
    virtual bool fileExists (const std::string& file) const
    {
      return db.fileExists(file);
    }


    /**
     * Open
     * @param file filename holding all data and keys
     * @param open_flags: can be regular UNIX file open flags such as: O_RDONLY,
     * O_RDWR, O_TRUNC
     * @param mode regular unix file mode
     *
     * @return 0 on success, -1 on failure with proper errno set
     *
     */
    virtual void open (const std::string& file, int open_flags, int mode)
    {
      int ret = db.open(file, open_flags, mode);

      if (ret != 0)
      {
	QDPIO::cerr << "Cannot open db file= " << file << std::endl;
	QDP_abort(1);
      }
    }


    virtual void close (void)
    {
      int ret = db.close();

      if (ret != 0)
      {
	QDPIO::cerr << "Error closing db file" << std::endl;
	QDP_abort(1);
      }
    }

    /**
     * Does this key exist in the store
     * @param key a key object
     * @return true if the answer is yes
     */
    bool exist(const K& key)
    {
      return db.exist(key);
    }


    /**
     * Insert a pair of data and key into the database
     * data is not ensemble, but a vector of complex.
     * @param key a key
     * @param data a user provided data
     *
     * @return 0 on successful write, -1 on failure with proper errno set
     */
    void insert (const K& key, const D& data)
    {
      int ret = db.insert(key, data);

      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error inserting into db" << std::endl;
	QDP_abort(1);
      }
    }


    /**
     * Insert a pair of data and key into the database in string format
     * @param key a key
     * @param data a user provided data
     *
     * @return 0 on successful write, -1 on failure with proper errno set
     */
    void insertBinary (const std::string& key, const std::string& data)
    {
      int ret = db.insertBinary(key, data);

      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error in db" << std::endl;
	QDP_abort(1);
      }
    }


    /**
     * Get data for a given key
     * @param key user supplied key
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int get (const K& key, D& data)
    {
      return db.get(key, data);
    }


    /**
     * Get data for a given key in binary form
     *
     * @param key user supplied key in string format
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int getBinary (const std::string& key, std::string& data)
    {
      return db.getBinary(key, data);
    }

    /**
     * Return all available keys to user
     * @param keys user suppled an empty vector which is populated
     * by keys after this call.
     */
    virtual void keys (std::vector<K>& keys_)
    {
      db.keys(keys_);
    }

    /**
     * Return all pairs of keys and data
     * @param keys user supplied empty vector to hold all keys
     * @param data user supplied empty vector to hold data
     * @return keys and data in the vectors having the same size
     */
    virtual void keysAndData (std::vector<K>& keys_, std::vector<D>& values_)
    {
      db.keysAndData(keys_, values_);
    }


    /**
     * Return all pairs of keys and data in binary string form
     * @param keys user supplied empty vector to hold all keys
     * @param data user supplied empty vector to hold data
     * @return keys and data in the vectors having the same size
     */
    virtual void binaryKeysAndData (std::vector<std::string> keys_,
				    std::vector<std::string> values_)
    {
      db.binaryKeysAndData(keys_, values_);
    }

    /**
     * Flush database in memory to disk
     */
    virtual void flush (void)
    {
      db.flush();
    }

    /**
     * Name of database associated with this Data store
     *
     * @return database name
     */
    virtual std::string storageName (void) const
    {
      return db.storageName();
    }


    /**
     * Insert user data into the  metadata database
     *
     * @param user_data user supplied data
     * @return returns 0 if success, else failure
     */
    virtual void insertUserdata (const std::string& user_data)
    {
      int ret = db.insertUserdata(user_data);

      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error inserting user data into db" << std::endl;
	QDP_abort(1);
      }
    }

    /**
     * Get user user data from the metadata database
     *
     * @param user_data user supplied buffer to store user data
     * @return returns 0 if success. Otherwise failure.
     */
    virtual void getUserdata (std::string& user_data)
    {
      int ret = db.getUserdata(user_data);

      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error getting user data from db" << std::endl;
	QDP_abort(1);
      }
    }


  private:
    FILEDB::ConfDataStoreDB<K, D> db;
  };


  //--------------------------------------------------------------------------------
  //!  DB Base class
  /*!
    This class is used for writing of user data (most usefully measurements)
    into a DB file with a key/value semantics.
  */
  template<typename K, typename D>
  class BinaryStoreDB
  {
  public:
    /**
     * Empty constructor for a DB
     */
    BinaryStoreDB ()
    {
      // Initialize default values
      if (Layout::primaryNode())
      {
	db.setCacheSize(db_cachesize);
//	db.enablePageMove();   // Disable!
	db.setMaxNumberConfigs(1);
      }
    }

    /*!
      Destroy the object
    */
    ~BinaryStoreDB()
    {
      close();
    }

    /**
     * How much data and keys should be kept in memory in bytes
     *
     * This should be called before the open is called
     * @param max_cache_size number of bytes of data and keys should be kept
     * in memory
     */
    virtual void setCacheSize (const unsigned int size)
    {
      if (Layout::primaryNode())
	db.setCacheSize(size);
    }


    /**
     * How much data and keys should be kept in memory in megabytes
     *
     * This should be called before the open is called
     * @param max_cache_size number of bytes of data and keys should be kept
     * in memory
     */
    virtual void setCacheSizeMB (const unsigned int size)
    {
      if (Layout::primaryNode())
	db.setCacheSizeMB(size);
    }


    /**
     * Page size used when a new data based is created
     * This only effects a new database
     *
     * @param pagesize the pagesize used in hash database. This value
     * should be power of 2. The minimum value is 512 and maximum value
     * is 262144
     *
     */
    virtual void setPageSize (const unsigned int size)
    {
      if (Layout::primaryNode())
	db.setPageSize(size);
    }


    /**
     * Set initial number of buckets
     *
     * This should be called before the open is called
     *
     * @param num the number of buckets should be used
     *
     */
    virtual void setNumberBuckets (const unsigned int num)
    {
      if (Layout::primaryNode())
	db.setNumberBuckets(num);
    }


    /**
     * Set whether to move pages when close to save disk space
     *
     * This only effective on writable database
     *
     */
    virtual void enablePageMove (void)
    {
      if (Layout::primaryNode())
	db.enablePageMove();
    }

    virtual void disablePageMove (void)
    {
      if (Layout::primaryNode())
	db.disablePageMove();
    }

    /**
     * Set and get maximum user information length
     */
    virtual void setMaxUserInfoLen (unsigned int len)
    {
      if (Layout::primaryNode())
	db.setMaxUserInfoLen(len);
    }

    virtual unsigned int getMaxUserInfoLen (void) const
    {
      unsigned int ret = 0;
      if (Layout::primaryNode())
	ret = db.getMaxUserInfoLen();

      QDPInternal::broadcast(ret);
      return ret;
    }

    /**
     * Set and get maximum number of configurations
     */
    virtual void setMaxNumberConfigs (unsigned int num)
    {
      if (Layout::primaryNode())
	db.setMaxNumberConfigs(num);
    }

    virtual unsigned int getMaxNumberConfigs (void) const
    {
      unsigned int ret = 0;
      if (Layout::primaryNode())
	ret = db.getMaxNumberConfigs();

      QDPInternal::broadcast(ret);
      return ret;
    }


    /**
     * Check if a DB file exists before opening.
     */
    virtual bool fileExists (const std::string& file) const
    {
      bool ret = false;
      if (Layout::primaryNode())
	ret = db.fileExists(file);

      QDPInternal::broadcast(ret);
      return ret;
    }


    /**
     * Open
     * @param file filename holding all data and keys
     * @param open_flags: can be regular UNIX file open flags such as: O_RDONLY,
     * O_RDWR, O_TRUNC
     * @param mode regular unix file mode
     *
     * @return 0 on success, -1 on failure with proper errno set
     *
     */
    virtual void open (const std::string& file, int open_flags, int mode)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.open(file, open_flags, mode);

      QDPInternal::broadcast(ret);
      if (ret != 0)
      {
	QDPIO::cerr << "Cannot open db file= " << file << std::endl;
	QDP_abort(1);
      }
    }


    virtual void close (void)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.close();

      QDPInternal::broadcast(ret);
      if (ret != 0)
      {
	QDPIO::cerr << "Error closing db file" << std::endl;
	QDP_abort(1);
      }
    }

    /**
     * Does this key exist in the store
     * @param key a key object
     * @return true if the answer is yes
     */
    bool exist(const K& key)
    {
      bool ret = true;
      if (Layout::primaryNode())
	ret = db.exist(key);

      QDPInternal::broadcast(ret);

      return ret;
    }


    /**
     * Insert a pair of data and key into the database
     * data is not ensemble, but a vector of complex.
     * @param key a key
     * @param data a user provided data
     *
     * @return 0 on successful write, -1 on failure with proper errno set
     */
    void insert (const K& key, const D& data)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.insert(key, data);

      QDPInternal::broadcast(ret);
      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error inserting into db" << std::endl;
	QDP_abort(1);
      }
    }


    /**
     * Insert a pair of data and key into the database in string format
     * @param key a key
     * @param data a user provided data
     *
     * @return 0 on successful write, -1 on failure with proper errno set
     */
    void insertBinary (const std::string& key, const std::string& data)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.insertBinary(key, data);

      QDPInternal::broadcast(ret);
      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error in db" << std::endl;
	QDP_abort(1);
      }
    }


    /**
     * Get data for a given key
     * @param key user supplied key
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int get (const K& key, D& data)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.get(key, data);
      else
	notImplemented();

      QDPInternal::broadcast(ret);
      return ret;
    }


    /**
     * Get data for a given key in binary form
     *
     * @param key user supplied key in string format
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int getBinary (const std::string& key, std::string& data)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.getBinary(key, data);
      else
	notImplemented();

      QDPInternal::broadcast(ret);
      return ret;
    }

    /**
     * Return all available keys to user
     * @param keys user suppled an empty vector which is populated
     * by keys after this call.
     */
    virtual void keys (std::vector<K>& keys_)
    {
      if (Layout::primaryNode())
	db.keys(keys_);
      else
	notImplemented();
    }

    /**
     * Return all pairs of keys and data
     * @param keys user supplied empty vector to hold all keys
     * @param data user supplied empty vector to hold data
     * @return keys and data in the vectors having the same size
     */
    virtual void keysAndData (std::vector<K>& keys_, std::vector<D>& values_)
    {
      if (Layout::primaryNode())
	db.keysAndData(keys_, values_);
      else
	notImplemented();
    }


    /**
     * Return all pairs of keys and data in binary string form
     * @param keys user supplied empty vector to hold all keys
     * @param data user supplied empty vector to hold data
     * @return keys and data in the vectors having the same size
     */
    virtual void binaryKeysAndData (std::vector<std::string> keys_,
				    std::vector<std::string> values_)
    {
      if (Layout::primaryNode())
	db.binaryKeysAndData(keys_, values_);
      else
	notImplemented();
    }

    /**
     * Flush database in memory to disk
     */
    virtual void flush (void)
    {
      if (Layout::primaryNode())
	db.flush();
    }

    /**
     * Name of database associated with this Data store
     *
     * @return database name
     */
    virtual std::string storageName (void) const
    {
      std::string ret;
      if (Layout::primaryNode())
	ret = db.storageName();

      QDPInternal::broadcast_str(ret);
      return ret;
    }


    /**
     * Insert user data into the  metadata database
     *
     * @param user_data user supplied data
     * @return returns 0 if success, else failure
     */
    virtual void insertUserdata (const std::string& user_data)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.insertUserdata(user_data);

      QDPInternal::broadcast(ret);
      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error inserting user data into db" << std::endl;
	QDP_abort(1);
      }
    }

    /**
     * Get user user data from the metadata database
     *
     * @param user_data user supplied buffer to store user data
     * @return returns 0 if success. Otherwise failure.
     */
    virtual void getUserdata (std::string& user_data)
    {
      int ret = 0;
      if (Layout::primaryNode())
	ret = db.getUserdata(user_data);
      /*
      else
        notImplemented();
      */

      QDPInternal::broadcast(ret);
      if (ret != 0)
      {
	QDPIO::cerr << __func__ << ": error getting user data from db" << std::endl;
	QDP_abort(1);
      }

      QDPInternal::broadcast_str(user_data);
    }


  private:
    void notImplemented() const
    {
      QDP_error_exit("FILEDB read routines do not work (yet) in parallel - only single node");
    }

    FILEDB::ConfDataStoreDB<K,D> db;
  };


}  // namespace Chroma

#endif
