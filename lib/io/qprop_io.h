/*
 *  Routines associated with simple propagator IO
 */

/*
 *  First the simple propagator header
 */

#ifndef __qprop_io_h__
#define __qprop_io_h__

struct PropHead{
  Real kappa;
  int source_smearingparam;
  int source_type;		// S-wave (0) or P-wave  (1)
  int source_direction;         // S-wave (0);   P-wave x(0) y(1) z(2)
  int sink_smearingparam;
  int sink_type;
  int sink_direction;
};

/*
 *  Routines for reading and writing propagator
 */

void readQprop(char file[], LatticePropagator& quark_prop, PropHead& header);

void writeQprop(char file[], const LatticePropagator& quark_prop, 
		const PropHead& header);

/*
 *  Routines for reading and writing the QQQ correlators
 */

void writeBarcomp(const char file[], const multiNd<Complex>& barprop, 
		  const PropHead& head_1, const PropHead& head_2,
		  const PropHead& head_3,
		  const int j_decay);

void readBarcomp(multiNd<Complex>& barprop, 
		 PropHead& head_1, PropHead& head_2, 
		 PropHead& head_3,
		 const char file[],  
		 const int j_decay);

// Routines for reading/writing headers

void writePropHead(const PropHead header, BinaryWriter& prop_out);
void writePropHeadinNml(const PropHead header, NmlWriter& nml);
void readPropHead(PropHead& header, BinaryReader& prop_in);


#endif
