// -*- C++ -*-
/*! \file
 *  \brief A Map Object that works lazily from Disk
 */


#ifndef __map_obj_disk_h__
#define __map_obj_disk_h__

#include "chromabase.h"
#include "util/ferm/map_obj.h"

using namespace QDP;

namespace Chroma
{

  class MapObjectDiskParams {
  public:

    // Constructor
    MapObjectDiskParams(XMLReader& xml_in, const std::string path) {
      XMLReader paramtop(xml_in, path);
      read(paramtop, "fileName", filename);
    }
    // Copy
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
  public:
    

    //! Default constructor
    MapObjectDisk(const MapObjectDiskParams& p_) : p(p_) {}


    //! OpenRead mode (Inserts) 
    bool openRead(void) {  
      bool ret_val=true;

      if ( !writer.is_open() ) {
	if( !reader.is_open() )  { 

	  // Open the reader
	  reader.open(p.getFileName());


	  /* FIND LOCATION OF METADATA, SEEK THERE, READ METADATA */


	  ret_val = true;
	}
	else { 
	  QDPIO::cerr << "MapObjectDisk: Already Open in read mode" << endl;
	  ret_val = true; // Re-opening in same mode is not an error
	}
	  
      }
      else { 
	QDPIO::cerr << "MapObjectDisk: Already  open in write mode" << endl;
	reader.close();
	ret_val = false; // ReadOpening without close from write mode is an error
      }
      
      return ret_val;
    }

    bool closeRead(void) { 
      bool ret_val = true;
      if( reader.is_open() ) {
	reader.close()
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

	/* Take note of current position. 
	   write out metadata dump
	   Rewind.
	   Skip MAGIC_NUMBER
	   write position of metadata
	   Seek end 
	   EOF
	*/
	   

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
      return (src_map.find(key) == src_map.end()) ? false : true;
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
	reader.seek(rpos);
	read(reader, val);
      }
      else { 
	QDPIO::cerr << "MapObjectDisk: Cant lookup if not open in read mode" << endl;
	throw std::string("MapObjectDisk: Cant lookup if not open in read mode");
      }
    }
			
    //! The number of elements
    typename MapType_t::size_type size() const {return src_map.size();}

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

    //! Map type convenience
    /* NB: An interesting question is whether I should use the 
       pos_type from BinaryReader or from BinaryWriter. BinaryWriter
       is which writes and BinaryReader is that which reads. Is there
       a unified std::pos_type? */
    typedef std::map<K, QDP::BinaryReader::pos_type> MapType_t;

    //! Usual begin iterator
    typename MapType_t::const_iterator begin() const {return src_map.begin();}

    //! Usual end iterator
    typename MapType_t::const_iterator end() const {return src_map.end();}

    //! Map of objects
    mutable MapType_t src_map;

    const MapObjectDiskParams p;
    BinaryFileReader reader;
    BinaryFileWriter writer;

  };

} // namespace Chroma

#endif
