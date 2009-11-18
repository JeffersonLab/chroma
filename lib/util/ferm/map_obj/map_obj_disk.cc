#include "map_obj_disk.h"
#include <string>

namespace Chroma { 
  namespace MapObjDiskEnv {
    const std::string file_magic="XXXXChromaLazyDiskMapObjFileXXXX";
  };

  MapObjDiskEnv::file_typenum_t 
  peekMapObjectDiskTypeCode(const std::string& filename)
  {
    BinaryFileReader reader;
    reader.open(filename);
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
    
    reader.close();
    return type_code;

  }
}
    
