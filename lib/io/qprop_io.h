/*
 *  Routines associated with simple propagator IO
 */

/*
 *  First the simple propagator header
 */

#ifndef __qprop_io_h__
#define __qprop_io_h__

//! Source header read
void read(XMLReader& xml, const string& path, PropHead& header);

//! Source header writer
void write(XMLWriter& xml, const string& path, const PropHead& header);

/*
 *  Routines for reading and writing propagator
 */


void readQprop(const string& file, LatticePropagator& quark_prop, PropHead& header);

void writeQprop(const string& file, const LatticePropagator& quark_prop, 
		const PropHead& header);

/*
 *  Routines for reading and writing the QQQ correlators
 */

void writeBarcomp(const string& file, const multiNd<Complex>& barprop, 
		  const PropHead& head_1, const PropHead& head_2,
		  const PropHead& head_3,
		  const int j_decay);

void readBarcomp(multiNd<Complex>& barprop, 
		 PropHead& head_1, PropHead& head_2, 
		 PropHead& head_3,
		 const string& file,  
		 const int j_decay);

// Routines for reading/writing headers

void writePropHead(const PropHead header, BinaryWriter& prop_out);
void writePropHeadinNml(const PropHead header, NmlWriter& nml);
void readPropHead(PropHead& header, BinaryReader& prop_in);


#endif
