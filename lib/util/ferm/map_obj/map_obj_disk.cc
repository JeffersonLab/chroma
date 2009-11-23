#include "map_obj_disk.h"
#include <string>

namespace Chroma { 
  namespace MapObjDiskEnv {
    // Magic string at start of file.
    const std::string file_magic="XXXXChromaLazyDiskMapObjFileXXXX";

    // Type hashes have to be shorter than 256 chars for now. 
    const size_t maxFileTypeLength=256;
  };

  // Disk traits.
  const string MapObjTraitsNum<KeyPropColorVec_t,LatticeFermion>::type_string="KeyTKeyPropColorVec_tValTLatticeFermion";
  const string MapObjTraitsNum<KeyBlockProp_t,LatticeFermion>::type_string="KeyTKeyBlockProp_tValTLatticeFermion";
  const string MapObjTraitsNum<KeyGridProp_t,LatticeFermion>::type_string="KeyTKeyGridProp_tValTLatticeFermion";
  const string MapObjTraitsNum<int,EVPair<LatticeColorVector> >::type_string="KeyTintValTEVPairLatticeColorVector";
  const string MapObjTraitsNum<char,float >::type_string="KeyTcharValTfloat";
  


  std::string
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
    
    std::string type_string;
    read(reader, type_string, MapObjDiskEnv::maxFileTypeLength);
#ifdef DISK_OBJ_DEBUGGING
    QDPIO::cout << "Read File Type String. Code=" << type_string << ". Current Position: " << reader.currentPosition() << endl;
#endif
    
    reader.close();
    return type_string;

  }
}
    
