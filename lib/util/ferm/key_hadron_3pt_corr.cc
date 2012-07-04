/*! \file
 * \brief Hadron 3pt correlators
 */

#include "util/ferm/key_hadron_3pt_corr.h"

namespace Chroma
{
#if 0
  namespace
  {
    //! Error output
    std::ostream& operator<<(std::ostream& os, const Array<int>& d)
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
    std::ostream& operator<<(std::ostream& os, const Array<ComplexD>& d)
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
  }


  //----------------------------------------------------------------------------
  //! Used for error output
  std::ostream& operator<<(std::ostream& os, const KeyHadron3PtCorr_t& d)
  {
    os << "KeyHadron3PtCorr_t:"
       << " src_name= " << d.src_name
       << " src_smear= " << d.src_smear
       << " src_lorentz= " << d.src_lorentz
       << " src_spin= " << d.src_spin
       << " snk_name= " << d.snk_name
       << " snk_smear= " << d.snk_smear
       << " snk_lorentz= " << d.snk_lorentz
       << " snk_spin= " << d.snk_spin
       << " mass= " << d.mass
       << " gamma= " << d.gamma
       << " links= " << d.links
       << " p_i= " << d.pi_pf.p_i
       << " p_f= " << d.pi_pf.p_f
       << " quark= " << d.quark
       << " dt= " << d.dt
       << " num_vecs= " << d.num_vecs
       << " ensemble= " << d.ensemble
       << std::endl;

    return os;
  }
#endif


  //----------------------------------------------------------------------------
  // Reader
  void read(XMLReader& xml, const std::string& path, PiPf& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "p_f", param.p_f);
    read(paramtop, "p_i", param.p_i);
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const PiPf& param)
  {
    push(xml, path);

    write(xml, "p_f", param.p_f);
    write(xml, "p_i", param.p_i);
    write(xml, "q_sq", norm2(param.p_f - param.p_i));

    pop(xml);
  }

  //----------------------------------------------------------------------------
  // Read a key
  void read(XMLReader& xml, const std::string& path, KeyHadron3PtCorr_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "num_vecs", param.num_vecs);
    read(paramtop, "src_name", param.src_name);
    read(paramtop, "src_smear", param.src_smear);
    read(paramtop, "src_lorentz", param.src_lorentz);
    read(paramtop, "src_spin", param.src_spin);
    read(paramtop, "snk_name", param.snk_name);
    read(paramtop, "snk_smear", param.snk_smear);
    read(paramtop, "snk_lorentz", param.snk_lorentz);
    read(paramtop, "snk_spin", param.snk_spin);
    read(paramtop, "PiPf", param.pi_pf);
    read(paramtop, "gamma", param.gamma);
    read(paramtop, "links", param.links);
    read(paramtop, "quark", param.quark);
    read(paramtop, "dt", param.dt);
    read(paramtop, "mass", param.mass);
    read(paramtop, "ensemble", param.ensemble);
  }

  // KeyHadron3PtCorr writer
  void write(XMLWriter& xml, const std::string& path, const KeyHadron3PtCorr_t& param)
  {
    push(xml, path);

    write(xml, "num_vecs", param.num_vecs);
    write(xml, "src_name", param.src_name);
    write(xml, "src_smear", param.src_smear);
    write(xml, "src_lorentz", param.src_lorentz);
    write(xml, "src_spin", param.src_spin);
    write(xml, "snk_name", param.snk_name);
    write(xml, "snk_smear", param.snk_smear);
    write(xml, "snk_lorentz", param.snk_lorentz);
    write(xml, "snk_spin", param.snk_spin);
    write(xml, "PiPf", param.pi_pf);
    write(xml, "gamma", param.gamma);
    write(xml, "links", param.links);
    write(xml, "quark", param.quark);
    write(xml, "dt", param.dt);
    write(xml, "mass", param.mass);
    write(xml, "ensemble", param.ensemble);

    pop(xml);
  }


  //----------------------------------------------------------------------------
  //! KeyHadron3PtCorr reader
  void read(BinaryReader& bin, KeyHadron3PtCorr_t& param)
  {
    read(bin, param.src_name, 128);
    read(bin, param.src_smear, 128);
    read(bin, param.src_lorentz);
    read(bin, param.src_spin);
    read(bin, param.snk_name, 128);
    read(bin, param.snk_smear, 128);
    read(bin, param.snk_lorentz);
    read(bin, param.snk_spin);
    read(bin, param.pi_pf.p_i);
    read(bin, param.pi_pf.p_f);
    read(bin, param.gamma);
    read(bin, param.links);
    read(bin, param.quark);
    read(bin, param.dt);
    read(bin, param.num_vecs);
    read(bin, param.mass, 128);
    read(bin, param.ensemble, 1024);
  }

  //! Hadron3PtCorr write
  void write(BinaryWriter& bin, const KeyHadron3PtCorr_t& param)
  {
    write(bin, param.src_name);
    write(bin, param.src_smear);
    write(bin, param.src_lorentz);
    write(bin, param.src_spin);
    write(bin, param.snk_name);
    write(bin, param.snk_smear);
    write(bin, param.snk_lorentz);
    write(bin, param.snk_spin);
    write(bin, param.pi_pf.p_i);
    write(bin, param.pi_pf.p_f);
    write(bin, param.gamma);
    write(bin, param.links);
    write(bin, param.quark);
    write(bin, param.dt);
    write(bin, param.num_vecs);
    write(bin, param.mass);
    write(bin, param.ensemble);
  }

} // namespace Chroma


