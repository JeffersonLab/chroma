/*! \file
 * \brief Key for distillution propagator sources and solutions
 */

#include <vector>
#include "util/ferm/key_prop_distillution.h"

namespace Chroma 
{ 
  namespace
  {
    //----------------------------------------------------------------------------
    //! a < b
    /*! This definition follows that of string comparison */
    template<typename T>
    inline bool operator<(const std::vector<T>& a, const std::vector<T>& b)
    {
      bool ret = false;
      int  len = (a.size() < b.size()) ? a.size() : b.size();

      for(int i=0; i < len; ++i)
      {
	if (a[i] != b[i])
	  return (a[i] < b[i]) ? true : false;
      }
    
      return (a.size() == b.size()) ? false : (a.size() < b.size()) ? true : false;
    }


    //----------------------------------------------------------------------------
    // Beastly hack - convert a string into a number sequence in an array
    std::vector<int> stringToArrayInt(const std::string& s)
    {
      std::vector<int> d;

      for(int i=0; i < s.size(); ++i)
	d.push_back(s[i]);

      return d;
    }

  } // end anonymous namespace


  //----------------------------------------------------------------------------
  // Support for the keys of prop colo
  bool operator<(const KeyPropDist_t& a, const KeyPropDist_t& b)
  {
    std::vector<int> lgaa;
    lgaa.push_back(stringToArrayInt(thePropDistTypeMap::Instance().lookUpString(a.line_type)));
    lgaa.push_back(a.t_source);
    lgaa.push_back(a.t_slice);
    lgaa.push_back(a.dist_src);
    lgaa.push_back(a.spin_src);
    lgaa.push_back(a.spin_snk);
    lgaa.push_back(a.quark_line);
    lgaa.push_back(stringToArrayInt(a.mass));

    std::vector<int> lgbb;
    lgbb.push_back(stringToArrayInt(thePropDistTypeMap::Instance().lookUpString(b.line_type)));
    lgbb.push_back(b.t_source);
    lgbb.push_back(b.t_slice);
    lgbb.push_back(b.dist_src);
    lgbb.push_back(b.spin_src);
    lgbb.push_back(b.spin_snk);
    lgbb.push_back(b.quark_line);
    lgbb.push_back(stringToArrayInt(b.mass));

    return (lgaa < lgbb);
  }



  //----------------------------------------------------------------------------
  // KeyPropDist read
  void read(BinaryReader& bin, KeyPropDist_t& param)
  {
    read(bin, param.line_type);
    read(bin, param.t_source);
    read(bin, param.t_slice);
    read(bin, param.dist_src);
    read(bin, param.spin_src);
    read(bin, param.spin_snk);
    read(bin, param.quark_line);
    read(bin, param.mass);
  }

  // KeyPropDist write
  void write(BinaryWriter& bin, const KeyPropDist_t& param)
  {
    write(bin, param.line_type);
    write(bin, param.t_source);
    write(bin, param.t_slice);
    write(bin, param.dist_src);
    write(bin, param.spin_src);
    write(bin, param.spin_snk);
    write(bin, param.quark_line);
    write(bin, param.mass);
  }

  //! KeyPropDist reader
  void read(XMLReader& xml, const std::string& path, KeyPropDist_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "line_type", param.line_type);
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "dist_src", param.dist_src);
    read(paramtop, "spin_src", param.spin_src);
    read(paramtop, "spin_snk", param.spin_snk);
    read(paramtop, "quark_line", param.quark_line);
    read(paramtop, "mass", param.mass);
  }

  // KeyPropDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropDist_t& param)
  {
    push(xml, path);

    write(xml, "line_type", param.line_type);
    write(xml, "t_source", param.t_source);
    write(xml, "t_slice", param.t_slice);
    write(xml, "dist_src", param.dist_src);
    write(xml, "spin_src", param.spin_src);
    write(xml, "spin_snk", param.spin_snk);
    write(xml, "quark_line", param.quark_line);
    write(xml, "mass", param.mass);

    pop(xml);
  }

} // namespace Chroma
