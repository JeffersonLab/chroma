/*
 *  $Id: qprop_io.cc,v 1.2 2003-04-17 20:09:38 dgr Exp $
 *
 *  These are a few simple I/O routines that we can use until QIO makes its appearance
 *  I have tried to include a simple header by means of a structure.
 */

#include "chromabase.h"
#include "common_io.h"
#include "io/qprop_io.h"
#include "primitives.h"
#include "qdp_util.h"

#include <string>
using std::string;

QDP_BEGIN_NAMESPACE(QDP);

/*
 *  A simple binary writer for the propagators
 *
 *  ARGUMENTS
 *  file: filename
 *  quark_prop: a quark propagator
 *  header: PropHead, defined in propagator.h
 */

void writeQprop(char file[], const LatticePropagator& quark_prop, 
		const PropHead& header){

  /*
   *  Now the local variables
   */
  
  BinaryWriter prop_out(file);

  writePropHead(header, prop_out); // Write out the header

  write(prop_out, quark_prop);	// The data
  
}

void readQprop(char file[], LatticePropagator& quark_prop, PropHead& header){

  /*
   *  Now the local variables
   */
  
  BinaryReader prop_in(file);
  
  readPropHead(header, prop_in); // Read the header

  read(prop_in, quark_prop);	// Now the data

}


void writeBarcomp(const char file[], const multiNd<Complex>& barprop, 
		  const PropHead& head_1, const PropHead& head_2,
		  const PropHead& head_3,
		  const int j_decay){

  /*
   * Local variables
   */
  
  int length = Layout::lattSize()[j_decay]; // Temporal extent of lattice

  multi1d<int> ranks(7);


  BinaryWriter qqq_out(file);	// Open the output file

  writePropHead(head_1, qqq_out);
  writePropHead(head_2, qqq_out);
  writePropHead(head_3, qqq_out); // Write ou the headers

  
  
  for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_3
    for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_2
      for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // sf_1
	for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_3
	  for(ranks[4]=0; ranks[4] < Ns; ++ranks[4])   // si_2
	    for(ranks[5]=0; ranks[5] < Ns; ++ranks[5]) // si_1
	      for(ranks[6] = 0; ranks[6] < length; ++ranks[6])
		{
		  write(qqq_out, barprop[ranks]);
		}
}


/*
 *  ~writePropHead writes a propagator header to an open binary file
 */

void writePropHead(const PropHead header, BinaryWriter& prop_out){

  /*
   *  The template for binary write performs a byte swap, and thus
   *  needs to know the type of the members of the structure.  Hence
   *  we cannot call a simple binary write on header.
   */

  write(prop_out, header.kappa);
  write(prop_out, header.source_smearingparam);
  write(prop_out, header.source_type);
  write(prop_out, header.source_direction);
  write(prop_out, header.sink_smearingparam);
  write(prop_out, header.sink_type);
  write(prop_out, header.sink_direction);
}

void writePropHeadinNml(const PropHead header, NmlWriter& nml) {
  Write(nml, header.kappa);
  Write(nml, header.source_smearingparam);
  Write(nml, header.source_type);
  Write(nml, header.source_direction);
  Write(nml, header.sink_smearingparam);
  Write(nml, header.sink_type);
  Write(nml, header.sink_direction);
}

void readPropHead(PropHead& header, BinaryReader& prop_in){

  /*
   *  The template for binary read performs a byte swap, and thus
   *  needs to know the type of the members of the structure.  Hence
   *  we cannot call a simple binary read on header.
   */

  read(prop_in, header.kappa);
  read(prop_in, header.source_smearingparam);
  read(prop_in, header.source_type);
  read(prop_in, header.source_direction);
  read(prop_in, header.sink_smearingparam);
  read(prop_in, header.sink_type);
  read(prop_in, header.sink_direction);
  
}

void readBarcomp(multiNd<Complex>& barprop, 
		 PropHead& head_1, PropHead& head_2, 
		 PropHead& head_3,
		 const char file[],  
		 const int j_decay){

  /*
   * Local variables
   */
  
  int length = Layout::lattSize()[j_decay]; // Temporal extent of lattice

  multi1d<int> ranks(7);


  BinaryReader qqq_in(file);	// Open the output file

  readPropHead(head_1, qqq_in); // Read the header
  readPropHead(head_2, qqq_in); // Read the header
  readPropHead(head_3, qqq_in); // Read the header

  
  
  for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_3
    for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_2
      for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // sf_1
	for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_3
	  for(ranks[4]=0; ranks[4] < Ns; ++ranks[4])   // si_2
	    for(ranks[5]=0; ranks[5] < Ns; ++ranks[5]) // si_1
	      for(ranks[6] = 0; ranks[6] < length; ++ranks[6])
		{
		  read(qqq_in, barprop[ranks]);
		}
}
