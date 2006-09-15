#include "actions/ferm/qprop/avp_inverter_interface.h"
using namespace QDP;

namespace Chroma { 
  
  namespace AVPSolverFunctions { 
    // Gauge Reader function - user supplied
      double gaugeReader(const void *OuterGauge,
			 void *env,
			 const int latt_coord[4],
			 int mu,
			 int row,
			 int col,
			 int reim)
      {
	/* Translate arg */
	multi1d<LatticeColorMatrix>& u = *(multi1d<LatticeColorMatrix>*)OuterGauge;
	
	// Get node and index
	multi1d<int> coord(Nd);
	coord = latt_coord;
	int node = Layout::nodeNumber(coord);
	int linear = Layout::linearSiteIndex(coord);
	
	if (node != Layout::nodeNumber()) {
	  
	  QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	  QDP_abort(1);
	}
	
	// Get the value
	// NOTE: it would be nice to use the "peek" functions, but they will
	// broadcast to all nodes the value since they are platform independent.
	// We don't want that, so we poke into the on-node data
	double val = (reim == 0) ? 
	  toDouble(u[mu].elem(linear).elem().elem(row,col).real()) : 
	  toDouble(u[mu].elem(linear).elem().elem(row,col).imag());
	
	return val;
      }
      
      
      // Fermion Reader function - user supplied
      double fermionReaderRHS(const void *OuterFermion,
			      void *env,
			      const int latt_coord[5],
			      int color,
			      int spin,
			      int reim)
      {
	/* Translate arg */
	multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)OuterFermion;
	int Ls1 = psi.size() - 1;
	
	// Get node and index
	int s = latt_coord[Nd];
	multi1d<int> coord(Nd);
	coord = latt_coord;
	int node = Layout::nodeNumber(coord);
	int linear = Layout::linearSiteIndex(coord);
	
	if (node != Layout::nodeNumber()) {
	  
	  QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	  QDP_abort(1);
	}
	
	// Get the value
	// NOTE: it would be nice to use the "peek" functions, but they will
	// broadcast to all nodes the value since they are platform independent.
	// We don't want that, so we poke into the on-node data
	double val = (reim == 0) ? 
	  double(psi[Ls1-s].elem(linear).elem(spin).elem(color).real()) : 
	  double(psi[Ls1-s].elem(linear).elem(spin).elem(color).imag());
	
	if (spin >= Ns/2)
	  val *= -1;
	
	return val;
      }
      
      
      // Fermion Reader function - user supplied
      double fermionReaderGuess(const void *OuterFermion,
				void *env,
				const int latt_coord[5],
				int color,
				int spin,
				int reim) 
      {
	/* Translate arg */
	multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)OuterFermion;
	int Ls1 = psi.size() - 1;
	
	// Get node and index
	int s = latt_coord[Nd];
	multi1d<int> coord(Nd);
	coord = latt_coord;
	int node = Layout::nodeNumber(coord);
	int linear = Layout::linearSiteIndex(coord);
	
	if (node != Layout::nodeNumber()) {
	  QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	  QDP_abort(1);
	}
	
	// Get the value
	// NOTE: it would be nice to use the "peek" functions, but they will
	// broadcast to all nodes the value since they are platform independent.
	// We don't want that, so we poke into the on-node data
	double val = (reim == 0) ? 
	  double(psi[Ls1-s].elem(linear).elem(spin).elem(color).real()) : 
	  double(psi[Ls1-s].elem(linear).elem(spin).elem(color).imag());
	
	if (spin >= Ns/2)
	  val *= -1;
	
	val *= -0.5;
	
	return val;
      }
      
      
      // Fermion Writer function - user supplied
      void fermionWriterSolver( void *OuterFermion, 
				void *env, 
				const int latt_coord[5],
				int color, 
				int spin,
				int reim,
				double val)
	
      {
	/* Translate arg */
	multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)OuterFermion;
	int Ls1 = psi.size() - 1;
	
	// Get node and index
	int s = latt_coord[Nd];
	multi1d<int> coord(Nd);
	coord = latt_coord;
	int node = Layout::nodeNumber(coord);
	int linear = Layout::linearSiteIndex(coord);
	
	if (node != Layout::nodeNumber()) {
	  QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	  QDP_abort(1);
	}
	
	// Rescale
	if (spin >= Ns/2)
	  val *= -1;
	
	val *= -2.0;
	
	// Set the value
	// NOTE: it would be nice to use the "peek" functions, but they will
	// broadcast to all nodes the value since they are platform independent.
	// We don't want that, so we poke into the on-node data
	if (reim == 0)
	  psi[Ls1-s].elem(linear).elem(spin).elem(color).real() = val;
	else
	  psi[Ls1-s].elem(linear).elem(spin).elem(color).imag() = val;
	
	return;
      }
      
      // Fermion Writer function - user supplied
      void fermionWriterOperator( void *OuterFermion, 
				  void *env, 
				  const int latt_coord[5],
				  int color, 
				  int spin,
				  int reim,
				  double val)
      {
	/* Translate arg */
	multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)OuterFermion;
	int Ls1 = psi.size() - 1;
	
	// Get node and index
	int s = latt_coord[Nd];
	multi1d<int> coord(Nd);
	coord = latt_coord;
	int node = Layout::nodeNumber(coord);
	int linear = Layout::linearSiteIndex(coord);
	
	if (node != Layout::nodeNumber()) {
	  QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	  QDP_abort(1);
	}
	
	// Rescale
	if (spin >= Ns/2)
	  val *= -1;
	
	val *= -0.5;
	
	// Set the value
	// NOTE: it would be nice to use the "peek" functions, but they will
	// broadcast to all nodes the value since they are platform independent.
	// We don't want that, so we poke into the on-node data
	if (reim == 0)
	  psi[Ls1-s].elem(linear).elem(spin).elem(color).real() = val;
	else
	  psi[Ls1-s].elem(linear).elem(spin).elem(color).imag() = val;
	
	return;
      }
  };

};
