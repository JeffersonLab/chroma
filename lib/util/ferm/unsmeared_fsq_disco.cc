/*! \file
 * \brief tr(F^2) disconnected blocks
 */

#include "util/ferm/unsmeared_fsq_disco.h"

#include <iostream>

namespace Chroma
{
  namespace
  {
   
    //! Error output
    std::ostream& operator<<(std::ostream& os, const multi1d<int>& d)
    {
      if (d.size() > 0)
      {
	os << d[0];

	for(int i=1; i < d.size(); ++i)
	  os << " " << d[i];
      }

      return os;
    }

    //! Error output
    std::ostream& operator<<(std::ostream& os, const std::vector<int>& d)
    {
      if (d.size() > 0)
      {
	os << d[0];

	for(int i=1; i < d.size(); ++i)
	  os << " " << d[i];
      }

      return os;
    }

    //! Error output
    std::ostream& operator<<(std::ostream& os, const std::vector< std::complex<double> >& d)
    {
      if (d.size() > 0)
      {
	for(int i=0; i < d.size(); ++i)
	{
	  if (i > 0)
	    os << "  ";

	  os << "(" << real(d[i]) << "," << imag(d[i]) << ")";
	}
      }

      return os;
    }

    //! Error output
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& d)
    {
      if (d.size() > 0)
      {
	os << d[0];

	for(int i=1; i < d.size(); ++i)
	  os << " " << d[i];
      }

      return os;
    }

    // Make a string
    std::string makeStr(const KeyFSqDiscoOperator_t& s1)
    {
      std::stringstream os;

      os << s1.smear << s1.left_lorentz << s1.right_lorentz << s1.disp_list << s1.mom;

      return os.str();
    }
  }



  //------------------------------------------------------------------------------
  // Comparison
  bool operator<(const KeyFSqDiscoOperator_t& s1, const KeyFSqDiscoOperator_t& s2)
  {
    return (makeStr(s1) < makeStr(s2));
  }
  
  //------------------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& os, const KeyFSqDiscoOperator_t& d)
  {
    os << "KeyFSqDiscoOperator_t:"
       << " smear= " << d.smear
       << " left_lorentz= " << d.left_lorentz
       << " right_lorentz= " << d.right_lorentz
       << " disp_list= " << d.disp_list
       << " mom= " << d.mom;

    return os;
  }

  
  //----------------------------------------------------------------------------
  // Read a key
  void read(XMLReader& xml, const std::string& path, KeyFSqDiscoOperator_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "smear", param.smear);
    read(paramtop, "left_lorentz", param.left_lorentz);
    read(paramtop, "right_lorentz", param.right_lorentz);
    read(paramtop, "disp_list", param.disp_list);
    read(paramtop, "mom", param.mom);
  }

  //! KeyHadron1PtCorr writer
  void write(XMLWriter& xml, const std::string& path, const KeyFSqDiscoOperator_t& param)
  {
    push(xml, path);
    
    write(xml, "smear", param.smear);
    write(xml, "left_lorentz", param.left_lorentz);
    write(xml, "right_lorentz", param.right_lorentz);
    write(xml, "disp_list", param.disp_list);
    write(xml, "mom", param.mom);

    pop(xml);
  }

  
  void write(XMLWriter& xml, const std::string& path, const ValFSqDiscoOperator_t& v)
  {
    push(xml, path);
    
    write(xml, "op", v.op);
    
    pop(xml); 
  }

  //! KeyOperator reader    
  void read(BinaryReader& bin, KeyFSqDiscoOperator_t& d)
  {
    readDesc(bin,d.smear);
    read(bin,d.left_lorentz);
    read(bin,d.right_lorentz);
    read(bin,d.disp_list);
    read(bin,d.mom);
  }
  
  //! KeyOperator writer
  void write(BinaryWriter& bin, const KeyFSqDiscoOperator_t& d)
  {
    writeDesc(bin,d.smear);
    write(bin,d.left_lorentz);
    write(bin,d.right_lorentz);
    write(bin,d.disp_list);
    write(bin,d.mom);
  }
  
  //! ValOperator reader    
  void read(BinaryReader& bin, ValFSqDiscoOperator_t& d)
  {
    read(bin,d.op);
  }
  
  //! ValOperator writer
  void write(BinaryWriter& bin, const ValFSqDiscoOperator_t& d)
  {
    write(bin,d.op);
  }

} // namespace Hadron


