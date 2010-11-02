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
    //! Concatenate two Array's
    template<typename T>
    inline std::vector<T> concat(const std::vector<T>& l, const std::vector<T>& r)
    {
      std::vector<int> nz(l.size() + r.size());
      int j = 0;
      for(int i=0; i < l.size(); ++i)
	nz[j++] = l[i];
  
      for(int i=0; i < r.size(); ++i)
	nz[j++] = r[i];

      return nz;
    }

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
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPropDist_t& param)
  {
    os << "KeyPropDist_t:";
    os << " prop_type= " << param.prop_type;
    os << " t_source= " << param.t_source;
    os << " t_slice= " << param.t_slice;
    os << " dist_src= " << param.dist_src;
    os << " spin_src= " << param.spin_src;
    os << " spin_snk= " << param.spin_snk;
    os << " quark_line= " << param.quark_line;
    os << " mass= " << param.mass;
    os << endl;

    return os;
  }

  //----------------------------------------------------------------------------
  // Support for the keys of prop colo
  bool operator<(const KeyPropDist_t& a, const KeyPropDist_t& b)
  {
    std::vector<int> lgaa;
    lgaa.push_back(a.t_source);
    lgaa.push_back(a.t_slice);
    lgaa.push_back(a.dist_src);
    lgaa.push_back(a.spin_src);
    lgaa.push_back(a.spin_snk);
    lgaa = concat(lgaa, stringToArrayInt(a.quark_line));
    lgaa = concat(lgaa, stringToArrayInt(a.prop_type));
    lgaa = concat(lgaa, stringToArrayInt(a.mass));

    std::vector<int> lgbb;
    lgbb.push_back(b.t_source);
    lgbb.push_back(b.t_slice);
    lgbb.push_back(b.dist_src);
    lgbb.push_back(b.spin_src);
    lgbb.push_back(b.spin_snk);
    lgbb = concat(lgbb, stringToArrayInt(b.quark_line));
    lgbb = concat(lgbb, stringToArrayInt(b.prop_type));
    lgbb = concat(lgbb, stringToArrayInt(b.mass));

    return (lgaa < lgbb);
  }



  //----------------------------------------------------------------------------
  // KeyPropDist read
  void read(BinaryReader& bin, KeyPropDist_t& param)
  {
    readDesc(bin, param.prop_type);
    read(bin, param.t_source);
    read(bin, param.t_slice);
    read(bin, param.dist_src);
    read(bin, param.spin_src);
    read(bin, param.spin_snk);
    readDesc(bin, param.quark_line);
    readDesc(bin, param.mass);
  }

  // KeyPropDist write
  void write(BinaryWriter& bin, const KeyPropDist_t& param)
  {
    writeDesc(bin, param.prop_type);
    write(bin, param.t_source);
    write(bin, param.t_slice);
    write(bin, param.dist_src);
    write(bin, param.spin_src);
    write(bin, param.spin_snk);
    writeDesc(bin, param.quark_line);
    writeDesc(bin, param.mass);
  }

  //! KeyPropDist reader
  void read(XMLReader& xml, const std::string& path, KeyPropDist_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "prop_type", param.prop_type);
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

    write(xml, "prop_type", param.prop_type);
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
