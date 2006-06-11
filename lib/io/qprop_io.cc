// $Id: qprop_io.cc,v 3.1 2006-06-11 06:30:33 edwards Exp $
/*! \file
 * \brief Routines associated with Chroma propagator IO
 */

#include "chromabase.h"
#include "io/param_io.h"
#include "io/qprop_io.h"

namespace Chroma 
{

  // Given a fermion action in string form, return the boundary
  /* HACK - THIS DEFINITELY NEEDS IMPROVEMENT */
  multi1d<int> getFermActBoundary(const string& fermact)
  {
    //
    // Initialize fermion action
    //
    std::istringstream  xml_s(fermact);
    XMLReader  fermacttop(xml_s);
    XMLReader  top(fermacttop, "/FermionAction");

    multi1d<int> boundary;

    try
    {
      if (top.count("FermionBC/boundary") != 0)
      {
	read(top, "FermionBC/boundary", boundary);
      }
      else if (top.count("boundary") != 0)
      {
	read(top, "boundary", boundary);
      }
      else
      {
	QDPIO::cerr << "Error: neither FermionBC group nor boundary found" << endl;
	throw string("Error: neither FermionBC group nor boundary found");
      }
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << "Error reading fermact: " << e << endl;
      throw e;
    }

    return boundary;
  }



  // Initialize header with default values
  PropSourceConst_t::PropSourceConst_t()
  {
    j_decay   = -1;
    t_source  = -1;
  }

  // Given a prop source xml in string form, return the t_srce
  multi1d<int> PropSourceConst_t::getTSrce()
  {
    //
    // Initialize source xml
    //
    std::istringstream  xml_s(source);
    XMLReader  sourcetop(xml_s);

    multi1d<int> t_srce;

    try
    {
      XMLReader  top(sourcetop, "/Source");

      read(top, "t_srce", t_srce);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << "Error reading source: " << e << endl;
      throw e;
    }

    return t_srce;
  }


  // Initialize header with default values
  PropSourceSmear_t::PropSourceSmear_t()
  {
    j_decay   = -1;
  }


  // Initialize header with default values
  PropSinkSmear_t::PropSinkSmear_t()
  {
    j_decay = -1;
  }


  // Initialize header with default values
  SeqSource_t::SeqSource_t()
  {
    j_decay   = -1;
    t_sink    = -1;
    sink_mom.resize(Nd-1);
    sink_mom = 0;
  }


  // Initialize header with default values
  ChromaProp_t::ChromaProp_t()
  {
    obsvP       = true;
    numRetries  = 1;
    // Create an document with an empty state info tag
  }

  // Initialize header with default values
  ChromaMultiProp_t::ChromaMultiProp_t()
  {
  }


  // Initialize header with default values
  QQQBarcomp_t::QQQBarcomp_t()
  {
    Dirac_basis = true;
    forward_props.resize(3);
  }


  // Initialize header with default values
  QQbarMescomp_t::QQbarMescomp_t()
  {
    Dirac_basis = true;
    forward_props.resize(3);
  }


  // Anonymous namespace
  namespace 
  {
    //! V5 Source header read 
    /*! This routine is SOLELY for backwards compatibility. It should go. */
    void readV5(XMLReader& xml, const string& path, PropSourceConst_t& header)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      read(paramtop, "source_type", header.source_type);

      std::string wave_state;
      read(paramtop, "wave_state", wave_state);
      if (wave_state != "S_WAVE")
      {
	throw std::string("version 5 only supports S_WAVE");
      }

      multi1d<int> t_srce;
      read(paramtop, "j_decay",  header.j_decay);
      read(paramtop, "t_source", t_srce);

      if (header.source_type == "SHELL_SOURCE")
      {
	XMLReader shelltop(paramtop, "ShellSource");

	std::string wvf_kind;
	Real wvf_param;
	int  wvfIntPar;
	{
	  XMLReader smeartop(shelltop, "SourceSmearingParam");
	    
	  read(smeartop, "wvf_kind", wvf_kind);
	  read(smeartop, "wvf_param", wvf_param);
	  read(smeartop, "wvfIntPar", wvfIntPar);
	}

	int laplace_power;
	read(shelltop, "laplace_power", laplace_power);
	if (laplace_power != 0)
	  throw std::string("only laplace_power=0 supported");

	Real link_smear_fact;
	int link_smear_num;
	read(shelltop, "link_smear_fact", link_smear_fact);
	read(shelltop, "link_smear_num", link_smear_num);

	int disp_length;
	int disp_dir;
	read(shelltop, "disp_length", disp_length);
	read(shelltop, "disp_dir", disp_dir);

	XMLReader xml_readback;
	{
	  XMLBufferWriter xml_tmp;
	  push(xml_tmp, "Param");
	  write(xml_tmp, "version", 6);
	  push(xml_tmp, "Source");
	  write(xml_tmp, "version", 1);
	  write(xml_tmp, "SourceType", header.source_type);

	  {
	    push(xml_tmp, "SmearingParam");
	    write(xml_tmp, "wvf_kind", wvf_kind);
	    write(xml_tmp, "wvf_param", wvf_param);
	    write(xml_tmp, "wvfIntPar", wvfIntPar);
	    write(xml_tmp, "no_smear_dir", header.j_decay);
	    pop(xml_tmp);
	  }

	  write(xml_tmp, "disp_length", disp_length);
	  write(xml_tmp, "disp_dir", disp_dir);

	  {
	    push(xml_tmp, "LinkSmearing");
	    write(xml_tmp, "LinkSmearingType", "APE_SMEAR");
	    write(xml_tmp, "link_smear_num", link_smear_num);
	    write(xml_tmp, "link_smear_fact", link_smear_fact);
	    write(xml_tmp, "no_smear_dir", header.j_decay);
	    pop(xml_tmp);
	  }

	  write(xml_tmp, "t_srce",  t_srce);
	  write(xml_tmp, "j_decay",  header.j_decay);

	  pop(xml_tmp);  // Source
	  pop(xml_tmp);  // Param

	  QDPIO::cout << "source_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	  xml_readback.open(xml_tmp);
	}

	// Recurse back in to re-read
	read(xml_readback, "/Param", header);
      }
      else if (header.source_type == "POINT_SOURCE")
      {
	XMLReader xml_readback;
	{
	  XMLBufferWriter xml_tmp;
	  push(xml_tmp, "Param");
	  write(xml_tmp, "version", 6);
	  push(xml_tmp, "Source");
	  write(xml_tmp, "version", 1);
	  write(xml_tmp, "SourceType", header.source_type);

	  write(xml_tmp, "t_srce",  t_srce);
	  write(xml_tmp, "j_decay",  header.j_decay);

	  pop(xml_tmp);  // Source
	  pop(xml_tmp);  // Param

	  QDPIO::cout << "source_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	  xml_readback.open(xml_tmp);
	}

	// Recurse back in to re-read
	read(xml_readback, "/Param", header);
      }
      else if (header.source_type == "WALL_SOURCE")
      {
	XMLReader xml_readback;
	{
	  XMLBufferWriter xml_tmp;
	  push(xml_tmp, "Param");
	  write(xml_tmp, "version", 6);
	  push(xml_tmp, "Source");
	  write(xml_tmp, "version", 1);
	  write(xml_tmp, "SourceType", header.source_type);

	  write(xml_tmp, "t_source",  t_srce[header.j_decay]);
	  write(xml_tmp, "j_decay",  header.j_decay);

	  pop(xml_tmp);  // Source
	  pop(xml_tmp);  // Param

	  QDPIO::cout << "source_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	  xml_readback.open(xml_tmp);
	}

	// Recurse back in to re-read
	read(xml_readback, "/Param", header);
      }
      else if (header.source_type == "RAND_Z2_WALL_SOURCE")
      {
	XMLReader xml_readback;
	{
	  XMLBufferWriter xml_tmp;
	  push(xml_tmp, "Param");
	  write(xml_tmp, "version", 6);
	  push(xml_tmp, "Source");
	  write(xml_tmp, "version", 1);
	  write(xml_tmp, "SourceType", header.source_type);

	  write(xml_tmp, "t_source",  t_srce[header.j_decay]);
	  write(xml_tmp, "j_decay",  header.j_decay);

	  pop(xml_tmp);  // Source
	  pop(xml_tmp);  // Param

	  QDPIO::cout << "source_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	  xml_readback.open(xml_tmp);
	}

	// Recurse back in to re-read
	read(xml_readback, "/Param", header);
      }
      else
      {
	std::ostringstream  os;
	os << "Unsupported source type= " << header.source_type 
	   << "  in file= " << __FILE__ 
	   << "  at line= " << __LINE__
	   << endl;
	throw os.str();
      }
    }


    //! V4 Sink header read 
    /*! This routine is SOLELY for backwards compatibility. It should go. */
    void readV4(XMLReader& xml, const string& path, PropSinkSmear_t& header)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      read(paramtop, "sink_type", header.sink_type);

      std::string wave_state;
      read(paramtop, "wave_state", wave_state);
      if (wave_state != "S_WAVE")
      {
	throw std::string("version 4 only supports S_WAVE");
      }

      read(paramtop, "j_decay", header.j_decay);

      if (header.sink_type == "SHELL_SINK")
      {
	XMLReader shelltop(paramtop, "ShellSink");

	std::string wvf_kind;
	Real wvf_param;
	int  wvfIntPar;
	{
	  XMLReader smeartop(shelltop, "SinkSmearingParam");
	    
	  read(smeartop, "wvf_kind", wvf_kind);
	  read(smeartop, "wvf_param", wvf_param);
	  read(smeartop, "wvfIntPar", wvfIntPar);
	}

	int laplace_power;
	read(shelltop, "laplace_power", laplace_power);
	if (laplace_power != 0)
	  throw std::string("only laplace_power=0 supported");

	Real link_smear_fact;
	int link_smear_num;
	read(shelltop, "link_smear_fact", link_smear_fact);
	read(shelltop, "link_smear_num", link_smear_num);

	int disp_length;
	int disp_dir;
	read(shelltop, "disp_length", disp_length);
	read(shelltop, "disp_dir", disp_dir);

	XMLReader xml_readback;
	{
	  XMLBufferWriter xml_tmp;
	  push(xml_tmp, "Param");
	  write(xml_tmp, "version", 6);
	  push(xml_tmp, "Sink");
	  write(xml_tmp, "version", 1);
	  write(xml_tmp, "SinkType", header.sink_type);

	  {
	    push(xml_tmp, "SmearingParam");
	    write(xml_tmp, "wvf_kind", wvf_kind);
	    write(xml_tmp, "wvf_param", wvf_param);
	    write(xml_tmp, "wvfIntPar", wvfIntPar);
	    write(xml_tmp, "no_smear_dir", header.j_decay);
	    pop(xml_tmp);
	  }

	  write(xml_tmp, "disp_length", disp_length);
	  write(xml_tmp, "disp_dir", disp_dir);

	  {
	    push(xml_tmp, "LinkSmearing");
	    write(xml_tmp, "LinkSmearingType", "APE_SMEAR");
	    write(xml_tmp, "link_smear_num", link_smear_num);
	    write(xml_tmp, "link_smear_fact", link_smear_fact);
	    write(xml_tmp, "no_smear_dir", header.j_decay);
	    pop(xml_tmp);
	  }

	  pop(xml_tmp);  // Sink
	  pop(xml_tmp);  // Param

	  QDPIO::cout << "sink_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	  xml_readback.open(xml_tmp);
	}

	// Recurse back in to re-read
	read(xml_readback, "/Param", header);
      }
      else if (header.sink_type == "POINT_SINK")
      {
	XMLReader xml_readback;
	{
	  XMLBufferWriter xml_tmp;
	  push(xml_tmp, "Param");
	  write(xml_tmp, "version", 6);
	  push(xml_tmp, "Sink");
	  write(xml_tmp, "version", 1);
	  write(xml_tmp, "SinkType", header.sink_type);

	  write(xml_tmp, "j_decay",  header.j_decay);

	  pop(xml_tmp);  // Sink
	  pop(xml_tmp);  // Param

	  QDPIO::cout << "sink_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	  xml_readback.open(xml_tmp);
	}

	// Recurse back in to re-read
	read(xml_readback, "/Param", header);
      }
      else if (header.sink_type == "WALL_SINK")
      {
	XMLReader xml_readback;
	{
	  XMLBufferWriter xml_tmp;
	  push(xml_tmp, "Param");
	  write(xml_tmp, "version", 6);
	  push(xml_tmp, "Sink");
	  write(xml_tmp, "version", 1);
	  write(xml_tmp, "SinkType", header.sink_type);

	  pop(xml_tmp);  // Sink
	  pop(xml_tmp);  // Param

	  QDPIO::cout << "sink_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	  xml_readback.open(xml_tmp);
	}

	// Recurse back in to re-read
	read(xml_readback, "/Param", header);
      }
      else
      {
	std::ostringstream  os;
	os << "Unsupported sink type= " << header.sink_type 
	   << "  in file= " << __FILE__ 
	   << "  at line= " << __LINE__
	   << endl;
	throw os.str();
      }
    }
  }



  // Source header read
  void read(XMLReader& xml, const string& path, PropSourceConst_t& header)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    try
    {
      switch (version) 
      {
      case 5:
	readV5(xml, path, header);
	break;

      case 6:
      {
	XMLReader xml_tmp(paramtop, "Source");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "SourceType", header.source_type);
	read(xml_tmp, "j_decay", header.j_decay);
	if (xml_tmp.count("t_source") != 0)
	{
	  read(xml_tmp, "t_source", header.t_source);
	}
	else if (xml_tmp.count("t_srce") != 0)
	{
	  multi1d<int> t_srce;
	  read(xml_tmp, "t_srce", t_srce);
	  if (t_srce.size() != Nd)
	  {
	    throw string("t_srce not expected size");
	  }
	  header.t_source = t_srce[header.j_decay];
	}
	else
	{
	  throw string("neither t_source nor t_srce found");
	}
	header.source = os.str();
      }
      break;

      default:
	/**************************************************************************/
	QDPIO::cerr << "PropSourceConst parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << "PropSourceConst: Error reading source: " << e << endl;
      QDP_abort(1);
    }
  }


  // Source header read
  void read(XMLReader& xml, const string& path, PropSourceSmear_t& header)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 6:
    {
      XMLReader xml_tmp(paramtop, "Source");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "SourceType", header.source_type);
      read(xml_tmp, "j_decay", header.j_decay);
      header.source = os.str();
    }
    break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "PropSourceSmear parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }
  }


  // Source header read
  void read(XMLReader& xml, const string& path, PropSinkSmear_t& header)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    try
    {
      switch (version) 
      {
      case 4:
	readV4(xml, path, header);
	break;

      case 5:
      {
	XMLReader xml_tmp(paramtop, "Sink");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "SinkType", header.sink_type); 
	read(xml_tmp, "j_decay", header.j_decay);
	header.sink = os.str();
      }
      break;

      default:
	/**************************************************************************/
	QDPIO::cerr << "PropSinkSmear parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << "PropSinkSmear: Error reading sink: " << e << endl;
      QDP_abort(1);
    }

  }

  //! SeqPropagator header reader
  void read(XMLReader& xml, const string& path, ChromaSeqProp_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
      // If you need previous versions, you'll have to go back to
      // the cvs source and dig for old versions. I'm not sure what
      // in production needs this support.

    case 3:
    {
      bool nonRelProp;
      read(paramtop, "nonRelSeqProp", nonRelProp); // new - is this prop non-relativistic
      if (nonRelProp)
	param.quarkSpinType = QUARK_SPIN_TYPE_UPPER;

      read(paramtop, "seq_src", param.seq_src);
    }
    break;

    case 4:
    {
      read(paramtop, "quarkSpinType", param.quarkSpinType);
      read(paramtop, "seq_src", param.seq_src);
    }
    break;

    default:
      QDPIO::cerr << "ChromaSeqProp parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "InvertParam", param.invParam);
    read(paramtop, "t_sink", param.t_sink);
    read(paramtop, "sink_mom", param.sink_mom);
  }


  //! SeqSource header reader
  void read(XMLReader& xml, const string& path, SeqSource_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
    {
      XMLReader xml_readback;
      {
	XMLBufferWriter xml_tmp;
	push(xml_tmp, "Param");
	write(xml_tmp, "version", 2);
	push(xml_tmp, "SeqSource");
	write(xml_tmp, "version", 1);
	std::string seq_src;
	read(paramtop, "seq_src", seq_src);
	write(xml_tmp, "SeqSourceType", seq_src);

	read(paramtop, "sink_mom", param.sink_mom);
	write(xml_tmp, "sink_mom",  param.sink_mom);

	read(paramtop, "t_sink", param.t_sink);
	write(xml_tmp, "t_sink",  param.t_sink);

	write(xml_tmp, "j_decay",  Nd-1);

	pop(xml_tmp);  // Source
	pop(xml_tmp);  // Param

	QDPIO::cout << "seqsrc_xml = XX" << xml_tmp.printCurrentContext() << "XX" << endl;

	xml_readback.open(xml_tmp);
      }

      // Recurse back in to re-read
      read(xml_readback, "/Param", param);
    }
    break;

    case 2:
    {
      XMLReader xml_tmp(paramtop, "SeqSource");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "SeqSourceType", param.seqsrc_type);
      read(xml_tmp, "t_sink", param.t_sink);
      read(xml_tmp, "sink_mom", param.sink_mom);
      read(xml_tmp, "j_decay", param.j_decay);
      param.seqsrc = os.str();
    }
    break;

    default:
      QDPIO::cerr << "SeqSource parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

  }


  // Forward propagator header read
  void read(XMLReader& xml, const string& path, ChromaProp_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    multi1d<int> boundary;

    param.quarkSpinType = QUARK_SPIN_TYPE_FULL;
    param.obsvP = true;
    param.numRetries = 1;

    QDPIO::cout << __func__ << ": numRetries = " << param.numRetries << endl;

    switch (version) 
    {
      /**************************************************************************/
    case 4:
    {
      // In V4 the fermion action specific stuff is within the <Param> tag and not
      // in a <FermionAction> tag beneath <Param>
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);

      XMLBufferWriter xml_out;
      push(xml_out,"FermionAction");
      xml_out << paramtop;
      pop(xml_out);

      param.fermact = xml_out.printCurrentContext();
    }
    break;

    /**************************************************************************/
    case 5:
    {
      // In this modified version of v4, the fermion action specific stuff
      // goes into a <FermionAction> tag beneath <Param>
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);

      XMLReader xml_tmp(paramtop, "FermionAction");

      XMLBufferWriter xml_out;
      push(xml_out,"FermionAction");
      xml_out << xml_tmp;
      write(xml_out, "boundary", boundary);
      pop(xml_out);

      param.fermact = xml_out.printCurrentContext();
    }
    break;

    /**************************************************************************/
    case 6:
    {
      bool nonRelProp;
      read(paramtop, "nonRelProp", nonRelProp); // new - is this prop non-relativistic
      if (nonRelProp)
	param.quarkSpinType = QUARK_SPIN_TYPE_UPPER;

      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);

      XMLReader xml_tmp(paramtop, "FermionAction");

      XMLBufferWriter xml_out;
      push(xml_out,"FermionAction");
      xml_out << xml_tmp;
      write(xml_out, "boundary", boundary);
      pop(xml_out);

      param.fermact = xml_out.printCurrentContext();
    }
    break;

    /**************************************************************************/
    case 7:
    {
      XMLReader xml_tmp(paramtop, "FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      param.fermact = os.str();

      read(paramtop, "InvertParam", param.invParam);

      bool nonRelProp;
      read(paramtop, "nonRelProp", nonRelProp); // new - is this prop non-relativistic
      if (nonRelProp)
	param.quarkSpinType = QUARK_SPIN_TYPE_UPPER;

      if (paramtop.count("obsvP") != 0)
      {
	read(paramtop, "obsvP", param.obsvP);
      }

      if (paramtop.count("boundary") != 0)
      {
	QDPIO::cerr << "ChromaProp: paranoia check - found a misplaced boundary" << endl; 
	QDP_abort(1);
      }
    }
    break;

    case 8:
    {
      XMLReader xml_tmp(paramtop, "FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      param.fermact = os.str();

      bool nonRelProp;
      read(paramtop, "nonRelProp", nonRelProp); // new - is this prop non-relativistic
      if (nonRelProp)
	param.quarkSpinType = QUARK_SPIN_TYPE_UPPER;

      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "obsvP", param.obsvP);

      if (paramtop.count("boundary") != 0)
      {
	QDPIO::cerr << "ChromaProp: paranoia check - found a misplaced boundary" << endl; 
	QDP_abort(1);
      }
    }
    break;

    case 9:
    {
      XMLReader xml_tmp(paramtop, "FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      param.fermact = os.str();

      read(paramtop, "quarkSpinType", param.quarkSpinType); // which quark spins to compute
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "obsvP", param.obsvP);

      if (paramtop.count("boundary") != 0)
      {
	QDPIO::cerr << "ChromaProp: paranoia check - found a misplaced boundary" << endl; 
	QDP_abort(1);
      }
    }
    break;

    case 10:
    {
      XMLReader xml_tmp(paramtop, "FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      param.fermact = os.str();

      read(paramtop, "quarkSpinType", param.quarkSpinType); // which quark spins to compute
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "obsvP", param.obsvP);
      read(paramtop, "numRetries", param.numRetries);
    
      if (paramtop.count("boundary") != 0)
      {
	QDPIO::cerr << "ChromaProp: paranoia check - found a misplaced boundary" << endl; 
	QDP_abort(1);
      }
    }
    break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "ChromaProp parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    QDPIO::cout << __func__ << ": exiting, param.numRetries = " << param.numRetries << endl;

  }


  // Forward propagator header read
  void read(XMLReader& xml, const string& path, ChromaMultiProp_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    param.quarkSpinType = QUARK_SPIN_TYPE_FULL;
    multi1d<int> boundary;

    switch (version) 
    {
      /**************************************************************************/
    case 5:  // Backward compatibility with non Multi ChromaProp_t
      /**************************************************************************/
    {
      read(paramtop, "MultiMasses", param.MultiMasses);
    
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);

      XMLReader xml_tmp(paramtop, "FermionAction");

      XMLBufferWriter xml_out;
      push(xml_out,"FermionAction");
      xml_out << xml_tmp;
      write(xml_out, "boundary", boundary);
      pop(xml_out);

      param.fermact = xml_out.printCurrentContext();
    }
    break;

    /**************************************************************************/
    case 6:  // Backward compatibility with non Multi ChromaProp_t
      /**************************************************************************/
    {
      read(paramtop, "MultiMasses", param.MultiMasses);

      bool nonRelProp;
      read(paramtop, "nonRelProp", nonRelProp); // new - is this prop non-relativistic
      if (nonRelProp)
	param.quarkSpinType = QUARK_SPIN_TYPE_UPPER;

      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);

      XMLReader xml_tmp(paramtop, "FermionAction");

      XMLBufferWriter xml_out;
      push(xml_out,"FermionAction");
      xml_out << xml_tmp;
      write(xml_out, "boundary", boundary);
      pop(xml_out);

      param.fermact = xml_out.printCurrentContext();
    }
    break;

    /**************************************************************************/
    case 7:
      /**************************************************************************/
    {
      read(paramtop, "MultiMasses", param.MultiMasses);

      XMLReader xml_tmp(paramtop, "FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      param.fermact = os.str();

      read(paramtop, "InvertParam", param.invParam);

      bool nonRelProp;
      read(paramtop, "nonRelProp", nonRelProp); // new - is this prop non-relativistic
      if (nonRelProp)
	param.quarkSpinType = QUARK_SPIN_TYPE_UPPER;

      if (paramtop.count("boundary") != 0)
      {
	QDPIO::cerr << "ChromaMultiProp: paranoia check - found a misplaced boundary" << endl; 
	QDP_abort(1);
      }
    }
    break;

    /* Compatibility with non MultiQprop */
    case 8:
    {
      read(paramtop, "MultiMasses", param.MultiMasses);

      XMLReader xml_tmp(paramtop, "FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      param.fermact = os.str();

      read(paramtop, "InvertParam", param.invParam);

      bool nonRelProp;
      read(paramtop, "nonRelProp", nonRelProp); // new - is this prop non-relativistic
      if (nonRelProp)
	param.quarkSpinType = QUARK_SPIN_TYPE_UPPER;
    }
    break;


    /* Compatibility with non MultiQprop */
    case 9:
    {
      read(paramtop, "MultiMasses", param.MultiMasses);

      XMLReader xml_tmp(paramtop, "FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      param.fermact = os.str();

      read(paramtop, "InvertParam", param.invParam);

      read(paramtop, "quarkSpinType", param.quarkSpinType); // quark spin components
    }
    break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "ChromaMultiProp parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    // Sanity check. No of residuals must be same as no of masses
    if ( param.MultiMasses.size() != param.invParam.RsdCG.size() ) { 
      QDPIO::cerr << "Number of Masses in param.MultiMasses (" 
		  << param.MultiMasses.size() 
		  << ") differs from no of RsdCGs in param.invParam.RsdCG (" 
		  << param.invParam.RsdCG.size() << ")" << endl;

      QDP_abort(1);
    }

  }


  //! SeqPropagator header reader
  void read(XMLReader& xml, const string& path, ForwardProp_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "PropSink", param.sink_header);
    read(paramtop, "ForwardProp", param.prop_header);
    read(paramtop, "PropSource", param.source_header);
  }


  //! SequentialSource header reader
  void read(XMLReader& xml, const string& path, SequentialSource_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "SeqSourceSinkSmear", param.sink_header);
    read(paramtop, "SeqSource", param.seqsource_header);
    read(paramtop, "ForwardProps", param.forward_props);
  }


  //! SequentialProp header reader
  void read(XMLReader& xml, const string& path, SequentialProp_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "SeqProp", param.seqprop_header);
    read(paramtop, "SeqSourceSinkSmear", param.sink_header);
    read(paramtop, "SeqSource", param.seqsource_header);
    read(paramtop, "ForwardProps", param.forward_props);
  }


  //! QQQBarcomp header reader
  void read(XMLReader& xml, const string& path, QQQBarcomp_t& param)
  {
    XMLReader paramtop(xml, path);

    int version = 1;
    if (paramtop.count("version") != 0)
      read(paramtop, "version", version);

    switch (version) 
    {
      /**************************************************************************/
    case 1:
      param.Dirac_basis = false;
      param.forward_props.resize(3);
      read(paramtop, "Propagator1", param.forward_props[0]);
      read(paramtop, "Propagator2", param.forward_props[1]);
      read(paramtop, "Propagator3", param.forward_props[2]);
      break;

      /**************************************************************************/
    case 2:
      read(paramtop, "Dirac_basis", param.Dirac_basis);
      read(paramtop, "ForwardProps", param.forward_props);
      break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "QQQBarcomp parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    if (param.forward_props.size() != 3)
    {
      QDPIO::cerr << "QQQBarcomp: unexpected number of forward_props = " 
		  << param.forward_props.size() << endl; 
      QDP_abort(1);
    }
  }


  //! QQbarMescomp header reader
  void read(XMLReader& xml, const string& path, QQbarMescomp_t& param)
  {
    XMLReader paramtop(xml, path);

    int version = 2;

    switch (version) 
    {
      /**************************************************************************/
    case 2:
      read(paramtop, "Dirac_basis", param.Dirac_basis);
      read(paramtop, "ForwardProps", param.forward_props);
      break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "QQbarMescomp parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    if (param.forward_props.size() != 2)
    {
      QDPIO::cerr << "QQbarMescomp: unexpected number of forward_props = " 
		  << param.forward_props.size() << endl; 
      QDP_abort(1);
    }
  }




  //---------------------------------------------------------------------------
  // Source header writer
  void write(XMLWriter& xml, const string& path, const PropSourceConst_t& header)
  {
    push(xml, path);

    int version = 6;
    write(xml, "version", version);
    xml << header.source;
    write(xml, "j_decay", header.j_decay);   // I think these two are duplicates of what
    write(xml, "t_source", header.t_source); // is in header.source

    pop(xml);
  }


  // Source header writer
  void write(XMLWriter& xml, const string& path, const PropSourceSmear_t& header)
  {
    push(xml, path);

    int version = 6;
    write(xml, "version", version);
    xml << header.source;
    write(xml, "j_decay", header.j_decay);

    pop(xml);
  }


  // Source header writer
  void write(XMLWriter& xml, const string& path, const PropSinkSmear_t& header)
  {
    push(xml, path);

    int version = 5;
    write(xml, "version", version);

    write(xml, "SinkType", header.sink_type);
    xml << header.sink;
    write(xml, "j_decay", header.j_decay);

    pop(xml);
  }


  // Write propagator inversion parameters
  void write(XMLWriter& xml, const string& path, const ChromaProp_t& header)
  {
    push(xml, path);

    int version = 10;
    write(xml, "version", version);
    write(xml, "quarkSpinType", header.quarkSpinType);
    write(xml, "obsvP", header.obsvP);           // new - measured 5D stuff
    write(xml, "numRetries", header.numRetries); 
    xml << header.fermact;
    write(xml, "InvertParam", header.invParam);

    pop(xml);
  }

  // Write propagator inversion parameters
  void write(XMLWriter& xml, const string& path, const ChromaMultiProp_t& header)
  {
    push(xml, path);

    int version = 8;
    write(xml, "version", version);
    write(xml, "quarkSpinType", header.quarkSpinType);
    write(xml, "MultiMasses", header.MultiMasses);
    xml << header.fermact;
    write(xml, "InvertParam", header.invParam);

    pop(xml);
  }


  //! SeqPropagator header writer
  void write(XMLWriter& xml, const string& path, const ChromaSeqProp_t& param)
  {
    push(xml, path);

    int version = 4;
    write(xml, "version", version);
    write(xml, "quarkSpinType", param.quarkSpinType);
    write(xml, "seq_src", param.seq_src);
    write(xml, "InvertParam", param.invParam);
    write(xml, "t_sink", param.t_sink);
    write(xml, "sink_mom", param.sink_mom);

    pop(xml);
  }


  //! SeqSource header writer
  void write(XMLWriter& xml, const string& path, const SeqSource_t& param)
  {
    push(xml, path);

    int version = 2;
    write(xml, "version", version);
    xml << param.seqsrc;

    pop(xml);
  }


  //! SeqPropagator header writer
  void write(XMLWriter& xml, const string& path, const ForwardProp_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "PropSinkSmear", param.sink_header);
    write(xml, "ForwardProp", param.prop_header);
    write(xml, "PropSource", param.source_header);

    if( path != "." )
      pop(xml);
  }


  //! SequentialSource header writer
  void write(XMLWriter& xml, const string& path, const SequentialSource_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "SeqSourceSinkSmear", param.sink_header);
    write(xml, "SeqSource", param.seqsource_header);
    write(xml, "ForwardProps", param.forward_props);

    if( path != "." )
      pop(xml);
  }


  //! SequentialProp header writer
  void write(XMLWriter& xml, const string& path, const SequentialProp_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "SeqProp", param.seqprop_header);
    write(xml, "SeqSourceSinkSmear", param.sink_header);
    write(xml, "SeqSource", param.seqsource_header);
    write(xml, "ForwardProps", param.forward_props);

    if( path != "." )
      pop(xml);
  }


  //! QQQBarcomp header writer
  void write(XMLWriter& xml, const string& path, const QQQBarcomp_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 2;
    write(xml, "version", version);
    write(xml, "Dirac_basis", param.Dirac_basis);
    write(xml, "ForwardProps", param.forward_props);

    if( path != "." )
      pop(xml);
  }


  //! QQbarMescomp header writer
  void write(XMLWriter& xml, const string& path, const QQbarMescomp_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 2;
    write(xml, "version", version);
    write(xml, "Dirac_basis", param.Dirac_basis);
    write(xml, "ForwardProps", param.forward_props);

    if( path != "." )
      pop(xml);
  }



  // Write a Chroma propagator
  /*
   * \param file_xml     file header ( Read )
   * \param record_xml   xml holding propagator info ( Read )
   * \param quark_prop   propagator ( Read )
   * \param file         path ( Read )
   * \param volfmt       either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
   * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
   */    
  void writeQprop(XMLBufferWriter& file_xml,
		  XMLBufferWriter& record_xml, const LatticePropagator& quark_prop,
		  const string& file, 
		  QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
  {
    QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
    write(to,record_xml,quark_prop);
    close(to);
  }


  // Write a Chroma propagator
  /*
   * \param file_xml     file header ( Read )
   * \param header       structure holding propagator info ( Read )
   * \param quark_prop   propagator ( Read )
   * \param file         path ( Read )
   * \param volfmt       either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
   * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
   */    
  void writeQprop(XMLBufferWriter& file_xml,
		  const ChromaProp_t& header, const LatticePropagator& quark_prop,
		  const string& file, 
		  QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
  {
    XMLBufferWriter record_xml;
    write(record_xml, "Propagator", header);
    writeQprop(file_xml, record_xml, quark_prop, file, volfmt, serpar);
  }




  // Read a Chroma propagator
  /*
   * \param file_xml     file header ( Write )
   * \param record_xml   xml holding propagator info ( Write )
   * \param quark_prop   propagator ( Write )
   * \param file         path ( Read )
   * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
   */    
  void readQprop(XMLReader& file_xml,
		 XMLReader& record_xml, LatticePropagator& quark_prop,
		 const string& file, 
		 QDP_serialparallel_t serpar)
  {
    QDPFileReader to(file_xml,file,serpar);
    read(to,record_xml,quark_prop);
    close(to);
  }

  // Read a Chroma propagator
  /*
   * \param file_xml     file header ( Write )
   * \param header       structure holding propagator info ( Write )
   * \param quark_prop   propagator ( Write )
   * \param file         path ( Read )
   * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
   */    
  void readQprop(XMLReader& file_xml,
		 ChromaProp_t& header, LatticePropagator& quark_prop,
		 const string& file, 
		 QDP_serialparallel_t serpar)
  {
    XMLReader record_xml;
    readQprop(file_xml, record_xml, quark_prop, file, serpar);
    read(record_xml, "/Propagator", header);  // extract header from xml
  }




  // Write a Chroma Fermion Field (eg prop_component)
  /*
   * \param file_xml     file header ( Read )
   * \param record_xml   xml holding propagator info ( Read )
   * \param fermion      fermion field( Read )
   * \param file         path ( Read )
   * \param volfmt       either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
   * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
   */    
  void writeFermion(XMLBufferWriter& file_xml,
		    XMLBufferWriter& record_xml, const LatticeFermion& fermion,
		    const string& file, 
		    QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
  {
    QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
    write(to,record_xml,fermion);
    close(to);
  }

  // Read a Chroma Fermion Field
  /*
   * \param file_xml     file header ( Write )
   * \param record_xml   xml holding propagator info ( Write )
   * \param fermion      The Fermion ( Write )
   * \param file         path ( Read )
   * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
   */    
  void readFermion(XMLReader& file_xml,
		   XMLReader& record_xml, 
		   LatticeFermion& fermion,
		   const string& file, 
		   QDP_serialparallel_t serpar)
  {
    QDPFileReader to(file_xml,file,serpar);
    read(to,record_xml,fermion);
    close(to);
  }

}  // end namespace Chroma

