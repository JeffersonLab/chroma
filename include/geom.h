// $Id: geom.h,v 1.2 2003-01-04 05:13:30 edwards Exp $
//
// QDP data parallel interface
//
// Geometry


//! Geometry namespace holding info on lattice
namespace Geometry
{
  //! Initializer for an isotropic geometry
  void init();

  //! Initializer for a anisotropic geometry
  void initAniso(int aniso_dir, float xi_0);

  //! Is anisotropy enabled?
  bool anisoP();

  //! Time direction
  int tDir();

  //! Anisotropy factor
  float xi_0();

};
