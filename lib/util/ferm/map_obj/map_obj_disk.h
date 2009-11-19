// -*- C++ -*-
/*! \file
 *  \brief A Map Object that works lazily from Disk
 */


#ifndef __map_obj_disk_h__
#define __map_obj_disk_h__

#include "chromabase.h"
#include "util/ferm/map_obj.h"
#include <string>
#include "util/ferm/map_obj/map_obj_disk_traits.h"

//#define DISK_OBJ_DEBUGGING 1
#undef DISK_OBJ_DEBUGGING

using namespace QDP;

namespace Chroma
{


  namespace MapObjDiskEnv { 
    typedef unsigned int file_version_t;
    typedef unsigned int file_typenum_t;

    extern const std::string file_magic;
  };


  //! A simple parameter class for MapObjectDisk types.
  class MapObjectDiskParams {
  public:
    MapObjectDiskParams() {};

    //! Constructor: From XML
    MapObjectDiskParams(XMLReader& xml_in, const std::string& path) {
      XMLReader paramtop(xml_in, path);
      read(paramtop, "FileName", filename);
    }

    //! Constructor: from filename
    MapObjectDiskParams(const std::string& inputFile): filename(inputFile) {} 

    //! Copy
    MapObjectDiskParams(const MapObjectDiskParams& p) : filename(p.filename) {}

    //! Destructor is automagic
    ~MapObjectDiskParams() {}

    //! const accessor 
    const std::string& getFileName() const { 
      return filename;
    }

    //! const mutator
    void setFileName(const std::string& _filename) { 
      filename = _filename;
    }
 
    void write(XMLWriter& xml_out, const std::string& path) const;
  private:
    std::string filename;
  };

#if 0
 void
 write(XMLWriter& xml_out, const std::string& path, const MapObjectDiskParams& p)
 {
   p.write(xml_out, path);
 }

 void
 read(XMLReader& xml_in, const std::string& path, MapObjectDiskParams& p) 
 {
   MapObjectDiskParams tmp(xml_in, path);
   p=tmp;
 }
#endif


  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObjectDisk : public MapObject<K,V>
  {
  public:
    

    //! Default constructor
    MapObjectDisk(const MapObjectDiskParams& p) : filename(p.getFileName()), file_version(1), state(INIT) {}

    MapObjectDisk(const std::string& filename_) : filename(filename_), file_version(1), state(INIT) {}

    /* --- STATE MACHINE FUNCTIONS: May Modify State -- */

    //! Open Write Mode
    void openWrite(void);

    //! Open Read mode (Lookups)
    void openRead(void);

    //! Open Update mode (Updates)
    void openUpdate(void);

    //! Finalizes object
    /*! Finalizes object from READ/WRITE states. 
     *  From any other state it throws and exception
     */
    ~MapObjectDisk();

    /* --- STATE AWARE FUNCTIONS: May require particular states -- */

    //! Insert: Requires WRITE state
    void insert(const K& key, const V& val);

    //! Lookup: Requires READ state
    void lookup(const K& key, V& val) const;

    //! Update: Requires UPDATE state
    void update(const K& key, const V& val);


    /* --- STATE AGNOSTIC METHODS: (work on in memory map structure -- */
    //! Exists?
    bool exist(const K& key) const;

    //! The number of elements
    unsigned int size() const {return static_cast<unsigned long>(src_map.size());}

    //! Dump keys
    std::vector<K> dump() const {
      std::vector<K> keys;
      if( reader.is_open() || writer.is_open() ) {
	
	typename MapType_t::const_iterator iter;
	for(iter  = src_map.begin();
	    iter != src_map.end();
	    ++iter) { 
	  keys.push_back(iter->first);
	}
      }
      return keys;

    }
    

  private:

    //! Type for the map
    typedef std::map<K, QDP::BinaryReader::pos_type> MapType_t;

    //! State machine states
    enum State { INIT, READ, WRITE, UPDATE };

    //! File related stuff. Unsigned int is as close to uint32 as I can get
    MapObjDiskEnv::file_version_t file_version;
    
    //! The state
    State state;
    
    //! Usual begin iterator
    typename MapType_t::const_iterator begin() const {return src_map.begin();}
    
    //! Usual end iterator
    typename MapType_t::const_iterator end() const {return src_map.end();}
    
    //! Map of objects
    mutable MapType_t src_map;
    
    //! The parameters
    const std::string filename;
    
    //! Reader and writer interfaces
    mutable BinaryFileReader reader;
    mutable BinaryFileWriter writer;
    
    // Internal Utility: Create/Skip past header
    void writeSkipHeader(void);
    
    //! Internal Utility: Read/Check header 
    BinaryReader::pos_type readCheckHeader(void);
    
    //! Internal Utility: Dump the map to disk
    void writeMapBinary(void);  
    
    //! Internal Utility: Read the map from disk
    void readMapBinary(const BinaryReader::pos_type& md_start);
    
    //! Internal Utility: Close File after write mode
    void closeWrite(void);
    
    //! Sink State for errors:
    void errorState(const std::string err) const {
      throw err;
    }


  };
  
  MapObjDiskEnv::file_typenum_t peekMapObjectDiskTypeCode(const std::string& filename);
  
  
  /* ****************** IMPLEMENTATIONS ********************** */
  
  /*! 
   * When called from INIT state, advances state machine to write mode.
   * Otherwise throws exception 
   */
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::openWrite(void) 
  {
    switch(state) { 
    case INIT: {
      QDPIO::cout << "MapObjectDisk: opening file " << filename
		  << " for writing" << endl;
      
      writer.open(filename);
          
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Writing file magic" << endl;
#endif
      // Write without newline or NULL -- Character by character
      const char* magic = MapObjDiskEnv::file_magic.c_str();
      writer.writeArray(magic,sizeof(char),MapObjDiskEnv::file_magic.length());
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Wrote magic. Current Position: " << writer.currentPosition() << endl;
#endif
      
      write(writer, (MapObjDiskEnv::file_version_t)file_version);
    
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Wrote Version. Current Position is: " << writer.currentPosition() << endl;
#endif
      
      MapObjDiskEnv::file_typenum_t type_code = MapObjTraitsNum<K,V>::filenum;
    
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Writing Type Code=" << type_code << endl;
#endif
      write(writer, type_code);
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Wrote Type Code. Current Position is: " << writer.currentPosition() << endl;
#endif
      
      BinaryReader::pos_type dummypos = static_cast<BinaryReader::pos_type>(writer.currentPosition());
    
#ifdef DISK_OBJ_DEBUGGING
      {
	QDPIO::cout << "Sanity Check 1" << endl; ;
	BinaryWriter::pos_type cur_pos = writer.currentPosition();
	if ( cur_pos != 
	     static_cast<BinaryReader::pos_type>(MapObjDiskEnv::file_magic.length())
	     +static_cast<BinaryReader::pos_type>(sizeof(MapObjDiskEnv::file_typenum_t))
	     +static_cast<BinaryReader::pos_type>(sizeof(MapObjDiskEnv::file_version_t)) ) {
	  
	  QDPIO::cout << "ERROR: Sanity Check 1 failed." << endl;
	  QDP_abort(1);
	}
      }
#endif
    
      /* Write a dummy link - make room for it */
      writer.writeArray((char *)&dummypos, sizeof(BinaryReader::pos_type), 1);
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Wrote dummy link: Current Position " << writer.currentPosition() << endl;
      {
	QDPIO::cout << "Sanity Check 2" << endl;
	BinaryWriter::pos_type cur_pos = writer.currentPosition();
	if ( cur_pos != 
	     static_cast<BinaryReader::pos_type>(MapObjDiskEnv::file_magic.length())
	     + static_cast<BinaryReader::pos_type>(sizeof(MapObjDiskEnv::file_version_t))
	     + static_cast<BinaryReader::pos_type>(sizeof(MapObjDiskEnv::file_typenum_t))
	     + static_cast<BinaryReader::pos_type>(sizeof(BinaryReader::pos_type)) ) {
	  QDPIO::cout << "Cur pos = " << cur_pos << endl;
	  QDPIO::cout << "Expected: " <<  MapObjDiskEnv::file_magic.length()+sizeof(MapObjDiskEnv::file_version_t) 
	    + sizeof(MapObjDiskEnv::file_typenum_t)
	    + sizeof(BinaryReader::pos_type) << endl;

	  QDPIO::cout << "ERROR: Sanity Check 2 failed." << endl;
	  QDP_abort(1);
	}
      }
#endif
      
      // Useful message
      QDPIO::cout << "MapObjectDisk: openWrite() Successful." << endl;
      
      // Advance state machine state
      state = WRITE;
      break;      
    }
    case UPDATE: {
      errorState("MapObjectDisk: Cannot openWrite() if in UPDATE mode");
      break;
    }
    
    case WRITE: {
      errorState("MapObjectDisk: Cannot openWrite() if already in write mode");
      break;
    }
    case READ:  {
      errorState("MapObjectDisk: Cannot openWrite() once in read mode");
      break;
    }
    
    default:
      errorState("MapOjectDisk: openWrite called from unknown state");
      break;
  }
  
    return;
}

 

  /*!
   * When called from READ state: does nothing. (Multiple clients can open read 
   * When called from INIT state: advances state machine to READ state.
   * When called from WRITE state: finalizes disk file and advances to READ 
   * When called from UPDATE state: closes Write file and resets to READ
   *   state.
   * otherwise it throws an exception 
   */
template<typename K, typename V>
  void
  MapObjectDisk<K,V>::openRead()
  {  
    switch (state) { 
    case UPDATE: {
      state = READ;
    }

    case READ:
      break ; // Multiple open reads are perfectly valid
              // they don't change data on disk or state
    case WRITE: {
      closeWrite(); // Complete the disk file
    }
      /*** DELIBERATE FALL THROUGH - to reopen in READ MODE ***/
    case INIT:
      {
	QDPIO::cout << "MapObjectDisk: opening file " << filename<< " for read access" << endl;
	
	// Open the reader
	reader.open(filename);
	
	QDPIO::cout << "MapObjectDisk: reading and checking header" << endl;
	BinaryReader::pos_type md_start = readCheckHeader();
	
	// Seek to metadata
	
	QDPIO::cout << "MapObjectDisk: reading key/fileposition data" << endl;
	
	/* Read the map in (metadata) */
	readMapBinary(md_start);
	
	/* And we are done */
	state = READ;
      }
      break;
    default:
      errorState("MapObjectDisk: openRead() called from unknown state");
      break;
    }
    
    return;
  }


  /*!
   * When called from READ state: opens Write file and sets update mode
   * When called from UPDATE state: does nothing
   * When called from INIT or WRITE states: throws exception
   */
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::openUpdate()
  {  
    switch (state) { 
    case UPDATE: // Do nothing. Deliberate Fallthrough
    case READ: {
      // Opens the writer -- will this blow away existing file
      state = UPDATE;
      break ; 
    }
    case WRITE: {
      errorState("MapObjectDisk: openUpdate() called from WRITE state.");
      break;
    }
    case INIT: {
      errorState("MapObjectDisk: openUpdate() called from INIT state.");
      break;
    }
    default:
      errorState("MapObjectDisk: openRead() called from unknown state");
      break;
    }
    
    return;
  }



  
  //! Destructor
  template<typename K, typename V>
  MapObjectDisk<K,V>::~MapObjectDisk() 
  {
    switch(state) { 
    case UPDATE: {
	writer.close();
	if( reader.is_open() ) { 
	  reader.close();
	}
	break;
    }

    case WRITE: {
	closeWrite(); // This finalizes files for us
	writer.close();
      }
      break;
    case READ:
      {
	reader.close();
	if( writer.is_open() ) { 
	  writer.close();
	}
      }
      break;
    default:
      // It is essentially useless to Destroy a Map in INIT state
      // It would finalize a file with NO DATA in it
      // Throw an exception.
      errorState("MapObjectDisk::Finalizing MAP Object Disk from INIT state");
      break;
    }
  }
  

  /*! 
   * Insert a value into the Map. Map has to be in WRITE state
   */
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::insert(const K& key, const V& val) 
  {
    switch (state)  { 
    case WRITE : {
      // Make note of current writer position
      BinaryWriter::pos_type pos=writer.currentPosition();
      
      // Turn writer pos into a reader pos
      BinaryReader::pos_type rpos = static_cast<BinaryReader::pos_type>(pos);
      
      // Add position to the map
      src_map.insert(std::make_pair(key,rpos));
      // Write object to disk
      write(writer, val);
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Wrote value to disk. Current Position: " << writer.currentPosition() << endl;
#endif
      break;
    }
    default:
      errorState("MapObjectDisk: insert attempted in non-write mode");
      break;
    }
  }

  /*! 
   * Update a value in the Map. Map has to be in UPDATE state
   */
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::update(const K& key, const V& val) 
  {
    switch (state)  { 
    case UPDATE : {
      //  Find key
      if (exist(key)){ 
	BinaryWriter::pos_type wpos = static_cast<BinaryWriter::pos_type>(src_map[key]);
#ifdef DISK_OBJ_DEBUGGING
       	QDPIO::cout << "Found key to update. Position is " << wpos << endl;
#endif
	
	
	writer.seek(wpos);
	
#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Sought write position. Current Position: " << writer.currentPosition() << endl;
#endif
	
	write(writer, val);

#ifdef DISK_OBJ_DEBUGGING	
	QDPIO::cout << "Wrote value to disk. Current Position: " << writer.currentPosition() << endl;
#endif
	writer.flush(); // Sync the file
      }
      else { 
	QDPIO::cerr << "MapObjectDisk: Cannot find key to update." << endl;
	errorState("MapObjectDisk: Cannot find key to update.");
      }
      break;
      
    }
    default:
      errorState("MapObjectDisk: insert attempted in non-write mode");
      break;
    }
  }


  /*! 
   * Lookup an item in the map.
   * Map has to be in write mode
   */
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::lookup(const K& key, V& val) const 
  { 
    switch(state) { 
    case UPDATE: // Deliberate fallthrough
    case READ: {
      
      if (! exist(key) ) {
	QDPIO::cerr << "MapObjectDisk: key not found" << std::endl;
	dump();
	errorState("MapObjectDisk: key not found");
      }
      
      // If key exists find file offset
      BinaryReader::pos_type rpos = src_map.find(key)->second;
      //      reader.rewind();
      reader.seek(rpos);
      read(reader, val);
      
      break;
    }
    default:
      errorState("MapObjectDisk: lookup() attempted when not in READ mode");
      break;
    }
  }
  
  // STATE AGNOSTIC FUNCTIONS (Work only on the in-memory-map state.
  template<typename K, typename V>
  bool 
  MapObjectDisk<K,V>::exist(const K& key) const 
  {
    return (src_map.find(key) == src_map.end()) ? false : true;
  }
  
  

  /***************** UTILITY ******************/


  //! Skip past header
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::writeSkipHeader(void) 
  { 
    if ( writer.is_open() ) { 
      writer.seek( MapObjDiskEnv::file_magic.length() + sizeof(MapObjDiskEnv::file_version_t) + sizeof(MapObjDiskEnv::file_typenum_t) );
    }
    else { 
      QDPIO::cerr << "Attempting writeSkipHeader, not in write mode" <<endl;
      QDP_abort(1);
    }
  }
  
  //! Check he header 
  template<typename K, typename V>
  BinaryReader::pos_type 
  MapObjectDisk<K,V>::readCheckHeader(void) {
    BinaryReader::pos_type md_position = 0;
    if( reader.is_open() ) {
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Rewinding File" << endl;
#endif
      
      reader.rewind();
      
      char* read_magic = new char[ MapObjDiskEnv::file_magic.length()+1 ];
      reader.readArray(read_magic, sizeof(char), MapObjDiskEnv::file_magic.length());
      read_magic[MapObjDiskEnv::file_magic.length()]='\0';
      
      // Check magic
      {
	std::string read_magic_str(read_magic);
	if (read_magic_str != MapObjDiskEnv::file_magic) { 
	  QDPIO::cerr << "Magic String Wrong: Expected: " << MapObjDiskEnv::file_magic << " but read: " << read_magic_str << endl;
	  QDP_abort(1);
	}
      }
      delete [] read_magic;

#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Read File Magic. Current Position: " << reader.currentPosition() << endl;
#endif
      
      MapObjDiskEnv::file_version_t read_version;
      read(reader, read_version);
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Read File Verion. Current Position: " << reader.currentPosition() << endl;
#endif
      
      // Check version
      QDPIO::cout << "MapObjectDisk: file has version: " << read_version << endl;
      
      MapObjDiskEnv::file_typenum_t type_code;
      read(reader, type_code);
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Read File Type Code. Code=" << type_code << ". Current Position: " << reader.currentPosition() << endl;
#endif
      
      // Read MD location
      reader.readArray((char *)&md_position, sizeof(BinaryReader::pos_type), 1);

#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Read MD Location. Current position: " << reader.currentPosition() << endl;
#endif
      
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Metadata starts at position: " << md_position << endl;
#endif
	
    }
    else { 
      QDPIO::cerr << "readCheckHeader needs reader mode to be opened. It is not" << endl;
      QDP_abort(1);
    }

    return md_position;
  }

  //! Dump the map
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::writeMapBinary(void)  
  {
    unsigned int map_size = src_map.size();
    write(writer, map_size);
#ifdef DISK_OBJ_DEBUGGING
    QDPIO::cout << "Wrote map size: " << map_size << " entries.  Current position : " << writer.currentPosition() << endl;
#endif
    
    typename MapType_t::const_iterator iter;
    for(iter  = src_map.begin();
	iter != src_map.end();
	++iter) { 
      
      K key = iter->first;
      BinaryReader::pos_type pos=iter->second;
      
      write(writer, key); 
      writer.writeArray((char *)&pos,sizeof(BinaryReader::pos_type),1);
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Wrote Key/Position pair:  Current Position: " << writer.currentPosition() << endl;
#endif
      
    }
  }
  
  //! read the map 
  // assume positioned at start of map data
  // Private utility function -- no one else should use.
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::readMapBinary(const BinaryReader::pos_type& md_start)
  {
    
    reader.seek(md_start);
#ifdef DISK_OBJ_DEBUGGING
    QDPIO::cout << "Sought start of metadata. Current position: " << reader.currentPosition() << endl;
#endif
    
    
    unsigned int num_records;
    
    read(reader, num_records);
#ifdef DISK_OBJ_DEBUGGING
    QDPIO::cout << "Read num of entries: " << num_records << " records. Current Position: " << reader.currentPosition() << endl;
#endif
    
    for(unsigned int i=0; i < num_records; i++) { 
      BinaryReader::pos_type rpos;
      K key;
      read(reader, key);
      
      reader.readArray( (char *)&rpos, sizeof(BinaryReader::pos_type),1);
      
#ifdef DISK_OBJ_DEBUGGING
      QDPIO::cout << "Read Key/Position pair. Current position: " << reader.currentPosition() << endl;
      
#endif
      // Add position to the map
      src_map.insert(std::make_pair(key,rpos));
    }
    
  }



  /*!
   * This is a utility function to sync the in memory offset map
   * with the one on the disk, and then close the file for writing.
   * Should not be called by user. Should only be called in WRITE state.
   */
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::closeWrite(void) 
{
    switch(state) { 
    case WRITE:
      {
	// Take note of current position
	BinaryReader::pos_type metadata_start =
	  static_cast<BinaryReader::pos_type>( writer.currentPosition() );
	
#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Beginning closeWrite: Metadata starts at position: "<<metadata_start << endl;
#endif
	// Dump metadata
	writeMapBinary();
	
	// Rewind and Skip header 
	writer.rewind();
#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Rewound file. Current Position: " << writer.currentPosition() << endl;
#endif
	writeSkipHeader();
#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Skipped Header. Current Position: " << writer.currentPosition() << endl;
#endif
	// write start position of metadata
	writer.writeArray((const char *)&metadata_start,sizeof(BinaryReader::pos_type),1);
	
#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Wrote link to metadata. Current Position: " << writer.currentPosition() << endl;
#endif
	
	// skip to end and close
	writer.seekEnd(0);
	writer.flush();
	
	QDPIO::cout << "MapObjectDisk: Closed file" << filename<< " for write access" <<  endl;
      }
      break;
    default:
      errorState("MapObjectDisk: closeFile() is an internal utility. Should only be called in WRITE mode");
      break;
    }
  }


} // namespace Chroma

#endif
