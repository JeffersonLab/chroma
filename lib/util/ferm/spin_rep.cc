/*! \file
 * \brief Sparse matrix representation of spin matrices
 */

#include "util/ferm/spin_rep.h"

namespace Chroma
{
  //----------------------------------------------------------------------------------
  // Construct a spin matrix in the DR rep
  SpinMatrix constructSpinDR(int gamma)
  {
    return Gamma(gamma) * SpinMatrix(Real(1));
  }

  //----------------------------------------------------------------------------------
  // Construct a spin matrix in the DP rep
  SpinMatrix constructSpinDP(int gamma)
  {
    return GammaDP(gamma) * SpinMatrix(Real(1));
  }


  //----------------------------------------------------------------------------
  //! MatrixSpinRep reader
  void read(XMLReader& xml, const std::string& path, MatrixSpinRep_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "left", param.left);
    read(paramtop, "right", param.right);
    read(paramtop, "Op", param.op);
  }

  //! MatrixSpinRep writer
  void write(XMLWriter& xml, const std::string& path, const MatrixSpinRep_t& param)
  {
    push(xml, path);

    write(xml, "left", param.left);
    write(xml, "right", param.right);
    write(xml, "Op", param.op);

    pop(xml);
  }

  //----------------------------------------------------------------------------
  //! MatrixSpinRep reader
  void read(BinaryReader& bin, MatrixSpinRep_t& param)
  {
    read(bin, param.left);
    read(bin, param.right);
    read(bin, param.op);
  }

  //! MatrixSpinRep writer
  void write(BinaryWriter& bin, const MatrixSpinRep_t& param)
  {
    write(bin, param.left);
    write(bin, param.right);
    write(bin, param.op);
  }


  //----------------------------------------------------------------------------------
  // Convert a generic spin matrix into a sparse spin representation
  std::vector<MatrixSpinRep_t> convertTwoQuarkSpin(const SpinMatrix& in)
  {
    std::vector<MatrixSpinRep_t> out;

    for(int sl = 0; sl < Ns; ++sl)
    {
      for(int sr = 0; sr < Ns; ++sr)
      {
	Complex ss = peekSpin(in, sl, sr);
	Real dd = localNorm2(ss);

	if (toBool(dd > Real(0.0)))
	{
	  MatrixSpinRep_t rep;
	  rep.left  = sl;
	  rep.right = sr;
	  rep.op    = ss;

	  out.push_back(rep);
	}
      }
    }

    return out;
  }

  //----------------------------------------------------------------------------------
  // Convert a DR gamma matrix indexed by gamma into a sparse spin representation
  std::vector<MatrixSpinRep_t> convertTwoQuarkSpinDR(int gamma)
  {
    return convertTwoQuarkSpin(constructSpinDR(gamma));
  }

  //----------------------------------------------------------------------------------
  // Convert a DP gamma matrix indexed by gamma into a sparse spin representation
  std::vector<MatrixSpinRep_t> convertTwoQuarkSpinDP(int gamma)
  {
    return convertTwoQuarkSpin(constructSpinDP(gamma));
  }


  //----------------------------------------------------------------------------------
  //! Fold in gamma_4 for source ops
  MatrixSpinRep_t foldSourceDP(const MatrixSpinRep_t& spin, bool src)
  {
    MatrixSpinRep_t obj(spin);

    int lfac = 1;
    int rfac = 1;

    if (src)
    {
      if (obj.left > 1)
	lfac = -1;
	
      if (obj.right > 1)
	rfac = -1;
    }

    obj.op *= Real(lfac * rfac);

    return obj;
  }

  //! Fold in gamma_4 for source ops
  std::vector<MatrixSpinRep_t> foldSourceDP(const std::vector<MatrixSpinRep_t>& spin, bool src)
  {
    std::vector<MatrixSpinRep_t> obj;

    for(int aks = 0; aks < spin.size(); ++aks)
    {
      obj.push_back(foldSourceDP(spin[aks], src));
    }

    return obj;
  }


  //----------------------------------------------------------------------------
  //! Rank3SpinRep reader
  void read(XMLReader& xml, const std::string& path, Rank3SpinRep_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "left", param.left);
    read(paramtop, "middle", param.middle);
    read(paramtop, "right", param.right);
    read(paramtop, "Op", param.op);
  }

  //! Rank3SpinRep writer
  void write(XMLWriter& xml, const std::string& path, const Rank3SpinRep_t& param)
  {
    push(xml, path);

    write(xml, "left", param.left);
    write(xml, "middle", param.middle);
    write(xml, "right", param.right);
    write(xml, "Op", param.op);

    pop(xml);
  }

  //----------------------------------------------------------------------------
  //! Rank3SpinRep reader
  void read(BinaryReader& bin, Rank3SpinRep_t& param)
  {
    read(bin, param.left);
    read(bin, param.middle);
    read(bin, param.right);
    read(bin, param.op);
  }

  //! Rank3SpinRep writer
  void write(BinaryWriter& bin, const Rank3SpinRep_t& param)
  {
    write(bin, param.left);
    write(bin, param.middle);
    write(bin, param.right);
    write(bin, param.op);
  }


  //----------------------------------------------------------------------------
  //! Convert a generic spin tensor into a sparse spin representation
  std::vector<Rank3SpinRep_t> convertThreeQuarkSpin(const Array3dO<Complex>& in)
  {
    std::vector<Rank3SpinRep_t> out;

    if ((in.size1() != Ns) || (in.size2() != Ns) || (in.size3() != Ns))
    {
      std::cerr << __func__ << ": input has unexpected size\n";
      exit(1);
    }

    for(int sl = 0; sl < Ns; ++sl)
    {
      for(int sm = 0; sm < Ns; ++sm)
      {
	for(int sr = 0; sr < Ns; ++sr)
	{
	  Complex ss = in(sl+1,sm+1,sr+1);
	  Real dd = localNorm2(ss);

	  // Avoid round-off funnies. Numbers should not be really small. Consider that zero.
	  if (toBool(dd > Real(1.0e-15)))
	  {
	    Rank3SpinRep_t rep;
	    rep.left   = sl;
	    rep.middle = sm;
	    rep.right  = sr;
	    rep.op     = ss;

	    out.push_back(rep);
	  }
	}
      }
    }

    return out;
  }


  //----------------------------------------------------------------------------------
  //! Fold in gamma_4 for source ops
  Rank3SpinRep_t foldSourceDP(const Rank3SpinRep_t& spin, bool src)
  {
    Rank3SpinRep_t obj(spin);

    int lfac = 1;
    int mfac = 1;
    int rfac = 1;

    if (src)
    {
      if (obj.left > 1)
	lfac = -1;
	
      if (obj.middle > 1)
	mfac = -1;
	
      if (obj.right > 1)
	rfac = -1;
    }

    obj.op *= Real(lfac * mfac * rfac);

    return obj;
  }

  //! Fold in gamma_4 for source ops
  std::vector<Rank3SpinRep_t> foldSourceDP(const std::vector<Rank3SpinRep_t>& spin, bool src)
  {
    std::vector<Rank3SpinRep_t> obj;

    for(int aks = 0; aks < spin.size(); ++aks)
    {
      obj.push_back(foldSourceDP(spin[aks], src));
    }

    return obj;
  }


}  // end namespace Chroma



  
