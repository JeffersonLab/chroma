// $Id: t_seqsource.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $
//! \file
//  \brief Test the sequential-source and resulting seqprop
//

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;

/*
 * Input 
 */
//! Parameters for running program
struct Param_t
{
  multi1d<int> nrow;
};

//! Propagators
struct Prop_t
{
  string           prop_file;  // The files is expected to be in SciDAC format!
  multi1d<string>  seqprop_files;  // The files is expected to be in SciDAC format!
};


//! Mega-structure of all input
struct Seqsource_input_t
{
  Param_t          param;
//  Cfg_t            cfg;
  Prop_t           prop;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "seqprop_files", input.seqprop_files);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "nrow", param.nrow);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Seqsource_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

//    // Read in the gauge configuration info
//    read(inputtop, "Cfg", input.cfg);

    // Read in the propagator(s) info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading prop data: " << e << endl;
    QDP_abort(1);
  }
}



int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  Seqsource_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/t_seqsource", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "t_seqsource");

  proginfo(xml_out);    // Print out basic program info

  //
  // Read the quark propagator and extract headers
  //
  XMLReader prop_file_xml, prop_record_xml;
  LatticePropagator quark_propagator;
  ChromaProp_t prop_header;
  PropSourceConst_t source_header;
  {
    QDPIO::cout << "Attempt to read forward propagator" << endl;
    readQprop(prop_file_xml, 
	      prop_record_xml, quark_propagator,
	      input.prop.prop_file, QDPIO_SERIAL);
   
    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try
    {
      read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
      read(prop_record_xml, "/Propagator/PropSource", source_header);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
      QDP_abort(1);
    }

    // Save propagator input
    write(xml_out, "Propagator_file_info", prop_file_xml);
    write(xml_out, "Propagator_record_info", prop_record_xml);
  }
  QDPIO::cout << "Forward propagator successfully read" << endl;

  // Derived from input prop
  int  j_decay = source_header.j_decay;
  multi1d<int> t_source = source_header.t_source;

  // Sanity check - write out the norm2 of the forward prop in the j_decay direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						 phases.getSet());

    push(xml_out, "Forward_prop_correlator");
    write(xml_out, "forward_prop_corr", forward_prop_corr);
    pop(xml_out);
  }



  XMLArrayWriter  xml_seq_src(xml_out, input.prop.seqprop_files.size());
  push(xml_seq_src, "Sequential_source");

  for (int seq_src_ctr = 0; seq_src_ctr < input.prop.seqprop_files.size(); ++seq_src_ctr) 
  {
    push(xml_seq_src);
    write(xml_seq_src, "seq_src_ctr", seq_src_ctr);

    // Read the sequential propagator
    // Read the quark propagator and extract headers
    LatticePropagator seq_quark_prop;
    SeqSource_t seqsource_header;
    {
      XMLReader seqprop_file_xml, seqprop_record_xml;
      readQprop(seqprop_file_xml, 
		seqprop_record_xml, seq_quark_prop,
		input.prop.seqprop_files[seq_src_ctr], QDPIO_SERIAL);

      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      // NEED SECURITY HERE - need a way to cross check props. Use the ID.
      try
      {
	read(seqprop_record_xml, "/SequentialProp/SeqSource", seqsource_header);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting seqprop header: " << e << endl;
	QDP_abort(1);
      }

      // Save seqprop input
      write(xml_seq_src, "SequentialProp_file_info", seqprop_file_xml);
      write(xml_seq_src, "SequentialProp_record_info", seqprop_record_xml);
    }
    QDPIO::cout << "Sequential propagator successfully read" << endl;

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);
      
      multi1d<Double> backward_prop_corr = sumMulti(localNorm2(seq_quark_prop), 
						    phases.getSet());
      
      push(xml_seq_src, "Backward_prop_correlator");
      write(xml_seq_src, "backward_prop_corr", backward_prop_corr);
      pop(xml_seq_src);
    }

    xml_out.flush();

    // Derived from input seqprop
    string seq_src = seqsource_header.seq_src;
    QDPIO::cout << "Seqsource name = " << seqsource_header.seq_src << endl;
    int           t_sink   = seqsource_header.t_sink;
    multi1d<int>  sink_mom = seqsource_header.sink_mom;

    write(xml_seq_src, "seq_src", seq_src);
    write(xml_seq_src, "t_source", t_source);
    write(xml_seq_src, "t_sink", t_sink);
    write(xml_seq_src, "sink_mom", sink_mom);
	
    // Now 

    int mom2_max = 0;
    for(int i=0; i < sink_mom.size(); ++i)
      mom2_max += sink_mom[i]*sink_mom[i];

    SftMom phases(mom2_max, sink_mom, false, j_decay);

    int G5 = Ns*Ns-1;
    int gamma_src = 0;
    int gamma_snk = 0;
    if (seq_src == "PION-A0_1")
    {
      gamma_src = G5;
      gamma_snk = 0;
    }
    else if (seq_src == "PION-A0_2")
    {
      gamma_src = G5;
      gamma_snk = 8;
    }
    else if (seq_src == "PION-RHO_X_1")
    {
      gamma_src = G5;
      gamma_snk = 1;
    }
    else if (seq_src == "PION-RHO_X_2")
    {
      gamma_src = G5;
      gamma_snk = 9;
    }
    else if (seq_src == "PION-RHO_Y_1")
    {
      gamma_src = G5;
      gamma_snk = 2;
    }
    else if (seq_src == "PION")
    {
      gamma_src = G5;
      gamma_snk = G5;
    }
    else if (seq_src == "PION-PION_2")
    {
      gamma_src = G5;
      gamma_snk = 7;
    }
    else
    {
      QDP_error_exit("unknown seq_src: %s", seq_src.c_str());
    }

    int gamma_insertion = G5 ^ gamma_src;

    push(xml_seq_src,"Corr_test");

    {
      Complex seq_src_corr = trace(peekSite(Gamma(G5)*adj(seq_quark_prop)*Gamma(G5)*Gamma(gamma_insertion), t_source));
      write(xml_seq_src,"gamma_insertion",gamma_insertion);
      write(xml_seq_src,"seq_src_corr",seq_src_corr);
    }

    {
      LatticePropagator anti_quark_prop =  Gamma(G5) * quark_propagator * Gamma(G5);
      LatticeComplex corr_fn = trace(adj(anti_quark_prop) * Gamma(gamma_snk) *
				     quark_propagator * Gamma(gamma_src));

      multi2d<DComplex> hsum;
      hsum = phases.sft(corr_fn);
//      for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num
//	multi1d<int> mom = phases.numToMom(sink_mom_num);

      Complex mesprop_0 = hsum[0][t_sink];
      write(xml_seq_src,"gamma_src",gamma_src);
      write(xml_seq_src,"gamma_snk",gamma_snk);
      write(xml_seq_src,"mesprop_0",mesprop_0);
    }

    pop(xml_seq_src);


    pop(xml_seq_src);   // elem
  } // end loop over sequential sources

  pop(xml_seq_src);  // Sequential_source

  pop(xml_out);     // t_seqsource

  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}

