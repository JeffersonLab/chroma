// $Id: geom.cc,v 1.3 2003-02-16 04:15:55 edwards Exp $
//
// QDP data parallel interface
//
// Geometry

#include "chroma.h"

namespace Geometry
{
  bool _aniso = false;;
  int  _t_dir = Nd-1;
  float _xi_0 = 1.0;

  //! Is anisotropy enabled?
  bool anisoP() {return _aniso;}

  //! Time direction
  int tDir() {return _t_dir;}

  //! Anisotropy factor
  float xi_0() {return _xi_0;}

  //! Initializer for geometry
  void init()
  {
    _xi_0 = 1.0;       // Anisotropy factor
    // By default, what is called the time direction is Nd-1
    // NOTE: nothing really depends on this except when aniso is turn on
    // The user can use any nrow direction for time
    _t_dir = Nd - 1;
  }


  //! Initializer for geometry
  void initAniso(int aniso_dir, float xx)
  {
    _aniso = true;   // No anisotropy by default
    _xi_0 = xx;       // Anisotropy factor
    _t_dir = aniso_dir;

    if (_t_dir < 0 || _xi_0 <= 0.0)
      QDP_error_exit("anisotropy values not set properly");
  }
};
