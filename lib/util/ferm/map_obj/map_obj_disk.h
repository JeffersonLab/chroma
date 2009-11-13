// -*- C++ -*-
/*! \file
 *  \brief A Map Object that works lazily from Disk
 */


#ifndef __map_obj_disk_h__
#define __map_obj_disk_h__

#include "chromabase.h"
#include "util/ferm/map_obj.h"
#include <string>

#define DISK_OBJ_DEBUGGING 1

using namespace QDP;

namespace Chroma
{

  
  class MapObjectDiskParams {
  public:

    //! Constructor: From XML
    MapObjectDiskParams(XMLReader& xml_in, const std::string path) {
      XMLReader paramtop(xml_in, path);
      read(paramtop, "fileName", filename);
    }

    //! Constructor: from filename
    MapObjectDiskParams(const std::string inputFile): filename(inputFile) {} 

    //! Copy
    MapObjectDiskParams(const MapObjectDiskParams& p) : filename(p.filename) {}

    // Destructor is automagic
    ~MapObjectDiskParams() {}

    const std::string& getFileName() const { 
      return filename;
    }

  private:
    std::string filename;
  };

  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObjectDisk : public MapObject<K,V>
  {
  private:
    typedef std::map<K, QDP::BinaryReader::pos_type> MapType_t;

  public:
    

    //! Default constructor
    MapObjectDisk(const MapObjectDiskParams& p_) : p(p_), file_magic(std::string("XXXXChromaLazyDiskMapObjFileXXXX")), file_version(1) {}


    //! OpenRead mode (Inserts) 
    bool openRead(void) {  
      bool ret_val=true;

      if ( !writer.is_open() ) {
	if( !reader.is_open() )  { 

	  // Open the reader
	  reader.open(p.getFileName());
	  BinaryReader::pos_type md_start = readCheckHeader();
	  
	  // Seek to metadata
	  reader.seek(md_start);

	  /* Read the map in (metadata) */
	  readMapBinary();

	  /* And we are done */
	  ret_val = true;
	}
	else { 
	  QDPIO::cerr << "MapObjectDisk: Already Open in read mode" << endl;
	  ret_val = true; // Re-opening in same mode is not an error
	}
	  
      }
      else { 
	QDPIO::cerr << "MapObjectDisk: Already  open in write mode" << endl;
	ret_val = false; // ReadOpening without close from write mode is an error
      }
      
      return ret_val;
    }

    bool closeRead(void) { 
      bool ret_val = true;
      if( reader.is_open() ) {
	reader.close();
	ret_val = true;
      }
      else {
	QDPIO::cerr << "MapObjectDisk: can\'t closeRead if not in read mode" << endl;
	ret_val = false; // 
      }

      return ret_val;
    }

    bool openWrite(void) {
      bool ret_val = true;
      if( !reader.is_open()  ) {
	if( !writer.is_open() ) { 
	  writer.open(p.getFileName());

	  /* Write a MAGIC NUMBER (filetype) and leave space for an
	     rpos type pointer */
#ifdef DISK_OBJ_DEBUGGING
	  QDPIO::cout << "Writing file magic" << endl;
#endif
	  // Write without newline or NULL -- Character by character
	  {
	    const char* magic = file_magic.c_str();
	    for(int tpos=0; tpos < file_magic.length(); tpos++)  {
	      write(writer, magic[tpos]);
	    }
	  }
	  QDPIO::cout << "Wrote magic. Position is now: " << writer.currentPosition() << endl;

#ifdef DISK_OBJ_DEBUGGING
	  QDPIO::cout << "Writing File version number" << endl;
#endif
	  write(writer, (unsigned long)file_version);

	  QDPIO::cout << "Wrote Version. Current position is: " << writer.currentPosition() << endl;

	  BinaryReader::pos_type dummypos = static_cast<BinaryReader::pos_type>(writer.currentPosition());


#ifdef DISK_OBJ_DEBUGGING
	  {
	    QDPIO::cout << "Sanity Check 1" << endl; ;
	    BinaryWriter::pos_type cur_pos = writer.currentPosition();
	    if ( cur_pos != 
		 file_magic.length()+sizeof(unsigned long) ) {
	      
	      QDPIO::cout << "ERROR: Sanity Check 1 failed." << endl;
	      QDP_abort(1);
	    }
	  }
#endif

	  /* Write a dummy link - make room for it */
	  writer.writeArray((char *)&dummypos, sizeof(BinaryReader::pos_type), 1);
	  QDPIO::cout << "Wrote dummy link: Current position " << writer.currentPosition() << endl;

#ifdef DISK_OBJ_DEBUGGING
	  {
	    QDPIO::cout << "Sanity Check 2" ;
	    BinaryWriter::pos_type cur_pos = writer.currentPosition();
	    if ( cur_pos != 
		 file_magic.length()+sizeof(unsigned long) + sizeof(BinaryReader::pos_type) ) {
	      QDPIO::cout << "Cur pos = " << cur_pos << endl;
	      QDPIO::cout << "Expected: " <<  file_magic.length()+sizeof(unsigned long) + sizeof(BinaryReader::pos_type) << endl;
	      QDPIO::cout << "ERROR: Sanity Check 2 failed." << endl;
	      QDP_abort(1);
	    }
	  }
#endif
	  /* We're done -- ready to start writing */
	  ret_val = true;
	}
	else { 
	  QDPIO::cerr << "MapObjectDisk: already open in write mode" << endl;
	  ret_val = true; // ReOpening from same mode is not really an error
	}
      }
      else { 
	QDPIO::cerr << "MapObjectDisk: can\'t OpenWrite if already open in read mode" << endl;
	ret_val = false;
      }
      return ret_val;
    }

    bool closeWrite(void) { 
      bool ret_val = true;
      if( writer.is_open() ) {
	
	/* Take note of current position. */
	BinaryReader::pos_type metadata_start =
	  static_cast<BinaryReader::pos_type>( writer.currentPosition() );

	/*
	 *  write out metadata dump
	 */

	/* Rewind and Skip header */
	writer.rewind();
	writeSkipHeader();

	/* write start position of metadata */
	write(writer, metadata_start);

	/* skip to end */
	writer.seekEnd(0);
	   
	/* Close */
	writer.close();
	
	ret_val = true;
      }
      else { 
	QDPIO::cerr << "Can\'t close write if not open in write mode" << endl;
	ret_val = false;
      }
      return ret_val;
    }

       

    //! Destructor
    ~MapObjectDisk() {
      if( reader.is_open()) closeRead();
      if( writer.is_open()) closeWrite();
    }

    //! Exists?
    bool exist(const K& key) const {
      if( reader.is_open() || writer.is_open() ) { 
	return (src_map.find(key) == src_map.end()) ? false : true;
      }
      else { 
	QDPIO::cerr << "exist called when neither in read/write mode. map may be uninitialized or incomplete" << endl;
      }
    }
			
    //! Insert
    void insert(const K& key, const V& val) {
      if(writer.is_open()) {
	
	// Make note of current writer position
	BinaryWriter::pos_type pos=writer.currentPosition();

	// Turn writer pos into a reader pos
	BinaryReader::pos_type rpos = static_cast<BinaryReader::pos_type>(pos);

	// Add position to the map
	src_map.insert(std::make_pair(key,rpos));

	// Write object to disk
	write(writer, val);
	
      }
      else { 
	QDPIO::cerr << "MapObjectDisk: Cant insert if not open in write mode" << endl;
	throw std::string("MapObjectDisk: Cant insert without being in write mode");
      }
    }
			
    //! Accessor
    void lookup(const K& key, V& val) const { 
      if( reader.is_open() ) {
	if (! exist(key) ) {
	  QDPIO::cerr << "MapObjectDisk: key not found" << std::endl;
	  // No generic key writer
	  //	QDPIO::cerr << "key= " << key << std::endl;
	  //	QDPIO::cerr << "All Keys:" << std::endl;
	  //	std::vector<K> all_keys = dump();
	  //	for(int i=0; i < all_keys.size(); ++i)
	  //	  QDPIO::cerr << all_keys[i];
	  
	  exit(1);
	}
	
	// If key exists find file offset
	BinaryReader::pos_type rpos = src_map.find(key)->second;
	reader.rewind();
	reader.seek(rpos);
	read(reader, val);
      }
      else { 
	QDPIO::cerr << "MapObjectDisk: Cant lookup if not open in read mode" << endl;
	throw std::string("MapObjectDisk: Cant lookup if not open in read mode");
      }
    }
			
    //! The number of elements
    unsigned long size() const {return static_cast<unsigned long>(src_map.size());}

    //! Dump keys
    std::vector<K> dump() const {
      std::vector<K> keys;
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter) { 
	keys.push_back(iter->first);
      }
      return keys;
    }
    

  private:

    std::string file_magic;
    unsigned long file_version;

    //! Map type convenience
    /* NB: An interesting question is whether I should use the 
       pos_type from BinaryReader or from BinaryWriter. BinaryWriter
       is which writes and BinaryReader is that which reads. Is there
       a unified std::pos_type? */


    //! Usual begin iterator
    typename MapType_t::const_iterator begin() const {return src_map.begin();}

    //! Usual end iterator
    typename MapType_t::const_iterator end() const {return src_map.end();}

    //! Map of objects
    mutable MapType_t src_map;

    const MapObjectDiskParams p;
    mutable BinaryFileReader reader;
    mutable BinaryFileWriter writer;

    //! Skip past header
    void writeSkipHeader(void) { 
      if ( writer.is_open() ) { 
	writer.seek( file_magic.length() + sizeof(unsigned long) );
      }
      else { 
	QDPIO::cerr << "Attempting writeSkipHeader, not in write mode" <<endl;
	QDP_abort(1);
      }
    }

    //! Check he header 
    BinaryReader::pos_type readCheckHeader(void) {
      BinaryReader::pos_type md_position = 0;
      if( reader.is_open() ) {
	reader.rewind();
	char* read_magic = new char[ file_magic.length()+1 ];
	for(int c=0; c < file_magic.length(); c++) { 
	  read(reader, read_magic[c]);
	}
	unsigned long read_version;
	read(reader, read_version);

	// Check magic
	read_magic[file_magic.length()]='\0';
	std::string read_magic_str(read_magic);
	if (read_magic_str != file_magic) { 
	  QDPIO::cerr << "Magic String Wrong: Expected: " << file_magic << " but read: " << read_magic_str << endl;
	  QDP_abort(1);
	}

	// Check version
	QDPIO::cout << "Read version#: " << read_version << endl;

	// No version dependent case statement for now


	// Read MD location
	char* md_pos = new char [ sizeof(BinaryReader::pos_type) ];
	reader.readArray(md_pos, sizeof(BinaryReader::pos_type), 1);
	md_position = static_cast< BinaryReader::pos_type >(*md_pos);
	delete [] md_pos;

	QDPIO::cout << "Metadata starts at position: " << md_position << endl;
	
      }
      else { 
	QDPIO::cerr << "readCheckHeader needs reader mode to be opened. It is not" << endl;
	QDP_abort(1);
      }
      return md_position;
    }

    //! Dump the map
    void writeMapBinary(void)  
    {
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter) { 
	write(writer, iter->first); 
	write(writer, iter->second);
      }
    }

    //! read the map 
    // assume positioned at start of map data
    // Private utility function -- no one else should use.
    void readMapBinary(void)
    {
      // Continue until we fail
      int reccount = 0;
      char* md_pos = new char [ sizeof(BinaryReader::pos_type ) ];

      while( !reader.fail() )  {
	BinaryReader::pos_type rpos;
	K key;
	read(reader, key);
	reader.readArray( md_pos, sizeof(BinaryReader::pos_type),1);
	rpos = static_cast<BinaryReader::pos_type>(*md_pos);
      

	// Add position to the map
	src_map.insert(std::make_pair(key,rpos));


	XMLBufferWriter buf_xml;
	push(buf_xml, "Record");
	write(buf_xml, "Key", key);
	pop(buf_xml);
	QDPIO::cout << "Read Record: "<< reccount
		    << " XML: " << buf_xml.str() << endl;


	reccount++;
      }
      
      QDPIO::cout << "Read " << reccount <<" records into map" << endl;
      delete [] md_pos;
    }

  };

} // namespace Chroma

#endif
