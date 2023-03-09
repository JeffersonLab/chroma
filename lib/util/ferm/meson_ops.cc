
#include "meson_ops.h"

namespace Chroma{
  namespace MesonOps{

    void read(XMLReader& xml, const std::string& path, Operators_t& op){
    XMLReader paramtop(xml, path);
    //QDPIO::cout<<"Reading Operators: "<<path<<std::endl ;
    //xml.print(std::cout);
    read(paramtop, "name"     ,op.name   );
    read(paramtop, "gamma"   ,op.gamma );
    }
    
    void write(XMLWriter& xml, const std::string& path, const Operators_t& op){
      push(xml, path);
      write(xml, "name"     ,op.name   );
      write(xml, "gamma"    ,op.gamma );
      pop(xml);
    }
    
    void read(XMLReader& xml,const std::string& path, State_t& s){
      XMLReader paramtop(xml, path);
      read(paramtop, "name"      ,s.name   );
      //read(paramtop, "key"       ,s.q      );
      read(paramtop, "db"        ,s.db     );
      read(paramtop, "wavefunc_file", s.wavefunc_file);
      //read(paramtop, "Operators" ,s.ops    );
    }
    
    void write(XMLWriter& xml, const std::string& path, const State_t& s){
      push(xml, path);
      write(xml, "name"      ,s.name   );
      //write(xml, "key"       ,s.q      );
      write(xml, "db"        ,s.db     );
      write(xml, "wavefunc_file", s.wavefunc_file);
      //write(xml, "Operators" ,s.ops    );
      pop(xml);
    }

  }

  namespace MesSpec{
    //! Key reader
    void read(XMLReader& xml, const std::string& path, Key& k){
      
      XMLReader paramtop(xml, path);
      
      read(paramtop, "parity",k.P );
      read(paramtop, "isospin",k.I );
      read(paramtop, "I3",k.I3 );
      read(paramtop, "spin",k.J );
      read(paramtop, "Jz",k.Jz );
      read(paramtop, "strangeness",k.S );
      read(paramtop, "charm",k.C );
      read(paramtop, "bottom",k.B );
      
    }
    
    //! Key writer
    void write(XMLWriter& xml, const std::string& path, const Key& k){
      
      push(xml, path);
      write(xml, "parity",k.P );
      write(xml, "isospin",k.I );
      write(xml, "I3",k.I3 );
      write(xml, "spin",k.J );
      write(xml, "Jz",k.Jz );
      write(xml, "strangeness",k.S );
      write(xml, "charm",k.C );
      write(xml, "bottom",k.B );
      pop(xml);
    }

  }// namespace MesSpec

}


