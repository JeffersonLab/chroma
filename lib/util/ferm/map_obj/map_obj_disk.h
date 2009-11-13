// -*- C++ -*-
/*! \file
 *  \brief A Map Object that works lazily from Disk
 */


#ifndef __map_obj_disk_h__
#define __map_obj_disk_h__

#include "chromabase.h"
#include "util/ferm/map_obj.h"
#include <string>

//#define DISK_OBJ_DEBUGGING 1
#undef DISK_OBJ_DEBUGGING

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

	  
	  QDPIO::cout << "MapObjectDisk: opening file " << p.getFileName() << " for read access" << endl;

	  // Open the reader
	  reader.open(p.getFileName());

	  QDPIO::cout << "MapObjectDisk: reading and checking header" << endl;
	  BinaryReader::pos_type md_start = readCheckHeader();
	  
	  // Seek to metadata

	  QDPIO::cout << "MapObjectDisk: reading key/fileposition data" << endl;

	  /* Read the map in (metadata) */
	  readMapBinary(md_start);

	  /* And we are done */
	  ret_val = true;
	  QDPIO::cout << "MapObjectDisk: openRead() successful" << endl;
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
	QDPIO::cout << "MapObjectDisk: closed file " << p.getFileName() <<" for read access" << endl;
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
	  QDPIO::cout << "MapObjectDisk: opening file " << p.getFileName() << " for writing" << endl;

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

#ifdef DISK_OBJ_DEBUGGING
	  QDPIO::cout << "Wrote magic. Current Position: " << writer.currentPosition() << endl;
#endif

	  write(writer, (file_version_t)file_version);
#ifdef DISK_OBJ_DEBUGGING
	  QDPIO::cout << "Wrote Version. Current Position is: " << writer.currentPosition() << endl;
#endif
	  BinaryReader::pos_type dummypos = static_cast<BinaryReader::pos_type>(writer.currentPosition());


#ifdef DISK_OBJ_DEBUGGING
	  {
	    QDPIO::cout << "Sanity Check 1" << endl; ;
	    BinaryWriter::pos_type cur_pos = writer.currentPosition();
	    if ( cur_pos != 
		 static_cast<BinaryReader::pos_type>(file_magic.length())
		 +static_cast<BinaryReader::pos_type>(sizeof(file_version_t)) ) {
	      
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
		 static_cast<BinaryReader::pos_type>(file_magic.length())
		 + static_cast<BinaryReader::pos_type>(sizeof(file_version_t))
		 + static_cast<BinaryReader::pos_type>(sizeof(BinaryReader::pos_type)) ) {
	      QDPIO::cout << "Cur pos = " << cur_pos << endl;
	      QDPIO::cout << "Expected: " <<  file_magic.length()+sizeof(file_version_t) + sizeof(BinaryReader::pos_type) << endl;
	      QDPIO::cout << "ERROR: Sanity Check 2 failed." << endl;
	      QDP_abort(1);
	    }
	  }
#endif
	  /* We're done -- ready to start writing */
	  QDPIO::cout << "MapObjectDisk: openWrite() Successful." << endl;
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

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Beginning closeWrite: Metadata starts at position: "<<metadata_start << endl;
#endif
	/*
	 *  write out metadata dump
	 */
	writeMapBinary();

	/* Rewind and Skip header */
	writer.rewind();
#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Rewound file. Current Position: " << writer.currentPosition() << endl;
#endif

	writeSkipHeader();

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Skipped Header. Current Position: " << writer.currentPosition() << endl;
#endif

	/* write start position of metadata */
	writer.writeArray((const char *)&metadata_start,sizeof(BinaryReader::pos_type),1);

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Wrote link to metadata. Current Position: " << writer.currentPosition() << endl;
#endif

	/* skip to end */
	writer.seekEnd(0);

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Sought end of file" << endl;
#endif 
	/* Close */
	writer.close();

	QDPIO::cout << "MapObjectDisk: Closed file" << p.getFileName() << " for write access" <<  endl;

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

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Wrote value to disk. Current Position: " << writer.currentPosition() << endl;
#endif
	
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
    typedef unsigned int file_version_t;

    std::string file_magic;
    file_version_t file_version;

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
	writer.seek( file_magic.length() + sizeof(file_version_t) );
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

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Rewinding File" << endl;
#endif

	reader.rewind();


	char* read_magic = new char[ file_magic.length()+1 ];
	reader.readArray(read_magic, sizeof(char), file_magic.length());
	read_magic[file_magic.length()]='\0';

	// Check magic
	{
	  std::string read_magic_str(read_magic);
	  if (read_magic_str != file_magic) { 
	    QDPIO::cerr << "Magic String Wrong: Expected: " << file_magic << " but read: " << read_magic_str << endl;
	    QDP_abort(1);
	  }
	}
	delete [] read_magic;

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Read File Verion. Current Position: " << reader.currentPosition() << endl;
#endif


	file_version_t read_version;
	read(reader, read_version);

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Read File Verion. Current Position: " << reader.currentPosition() << endl;
#endif


	// Check version
	QDPIO::cout << "MapObjectDisk: file has version: " << read_version << endl;

	// No version dependent case statement for now


	// Read MD location
	char* md_pos = new char [ sizeof(BinaryReader::pos_type) ];
	reader.readArray(md_pos, sizeof(BinaryReader::pos_type), 1);

#ifdef DISK_OBJ_DEBUGGING
	QDPIO::cout << "Read MD Location. Current position: " << reader.currentPosition() << endl;
#endif

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
    void readMapBinary(const BinaryReader::pos_type& md_start)
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

  };

} // namespace Chroma

#endif
