// $Id: geom.h,v 1.1 2004-01-08 03:11:47 edwards Exp $
//
// QDP data parallel interface
//
// Geometry

#error "NOT CURRENTLY USED."

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
