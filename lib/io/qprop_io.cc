// $Id: qprop_io.cc,v 1.30 2005-03-07 02:54:15 edwards Exp $
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
  PropSource_t::PropSource_t()
  {
    source_type = SRC_TYPE_POINT_SOURCE;
    wave_state  = WAVE_TYPE_S_WAVE;

    j_decay         = 0;
    direction       = 0;
    laplace_power   = 0;
    disp_length     = 0;
    disp_dir        = 0;
    link_smear_fact = 0;
    link_smear_num  = 0;
    t_source.resize(Nd);
    t_source        = 0;
    nrow            = Layout::lattSize();
  }

  // Initialize header with default values
  PropSink_t::PropSink_t()
  {
    sink_type     = SNK_TYPE_POINT_SINK;
    wave_state    = WAVE_TYPE_S_WAVE;

    direction       = 0;
    laplace_power   = 0;
    disp_length     = 0;
    disp_dir        = 0;
    link_smear_fact = 0;
    link_smear_num  = 0;
    nrow            = Layout::lattSize();
  }

  // Initialize header with default values
  PropSink_t::PropSink_t(const PropSource_t& source)
  {
    // Convert the source to a sink type
    switch(source.source_type)
    {
    case SRC_TYPE_POINT_SOURCE:
      sink_type = SNK_TYPE_POINT_SINK;
      break;
    case SRC_TYPE_WALL_SOURCE:
      sink_type = SNK_TYPE_WALL_SINK;
      break;
    case SRC_TYPE_SHELL_SOURCE:
      sink_type = SNK_TYPE_SHELL_SINK;
      break;
    case SRC_TYPE_BNDST_SOURCE:
      sink_type = SNK_TYPE_BNDST_SINK;
      break;
    case SRC_TYPE_POINT_AND_BNDST_SOURCE:
      sink_type = SNK_TYPE_POINT_AND_BNDST_SINK;
      break;
    case SRC_TYPE_SHELL_AND_BNDST_SOURCE:
      sink_type = SNK_TYPE_SHELL_AND_BNDST_SINK;
      break;
    case SRC_TYPE_POINT_AND_SHELL_AND_BNDST_SOURCE:
      sink_type = SNK_TYPE_POINT_AND_SHELL_AND_BNDST_SINK;
      break;
    }

    // description (same order):
    // wave_state: "S_WAVE", "P_WAVE", or "D_WAVE"
    // smearing params
    //   smearing width
    //   number of iteration for smearing
    //   number of iteration for smearing
    // decay direction
    // power of laplacian operator
    // displacement length
    // displacement direction: x(0),y(1),z(2)
    // sink_dir: direction of derivative at sink
    // smearing factor
    // number of smearing hits
    // Maximum number of blocking/smearing iterations
    // Blocking/smearing accuracy

    wave_state      = source.wave_state;
    sinkSmearParam  = source.sourceSmearParam;
    direction       = source.direction;
    laplace_power   = source.laplace_power;
    disp_length     = source.disp_length;
    disp_dir        = source.disp_dir;
    link_smear_fact = source.link_smear_fact;
    link_smear_num  = source.link_smear_num;
    nrow            = Layout::lattSize();
  }



  // Initialize header with default values
  ChromaProp_t::ChromaProp_t()
  {
    nrow        = Layout::lattSize();
    // Create an document with an empty state info tag
  }

  // Initialize header with default values
  ChromaMultiProp_t::ChromaMultiProp_t()
  {
    nrow        = Layout::lattSize();
  }


  // Initialize header with default values
  QQQBarcomp_t::QQQBarcomp_t()
  {
    Dirac_basis = true;
    forward_props.resize(3);
  }



  // Source header read
  void read(XMLReader& xml, const string& path, PropSource_t& header)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
      /**************************************************************************/
    case 5:
      /**************************************************************************/
      break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "PropSource parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "wave_state", header.wave_state);
    if (header.wave_state != WAVE_TYPE_S_WAVE)
    {
      read(paramtop, "direction",  header.direction);
    }

    read(paramtop, "source_type", header.source_type);
    if (header.source_type == SRC_TYPE_SHELL_SOURCE)
    {
      XMLReader shelltop(paramtop, "ShellSource");

      read(shelltop, "SourceSmearingParam", header.sourceSmearParam);
      read(shelltop, "laplace_power", header.laplace_power);
      read(shelltop, "link_smear_fact", header.link_smear_fact);
      read(shelltop, "link_smear_num", header.link_smear_num);
      read(shelltop, "disp_length", header.disp_length);
      read(shelltop, "disp_dir", header.disp_dir);
    }		

    read(paramtop, "j_decay",  header.j_decay);
    read(paramtop, "t_source", header.t_source);

    read(paramtop, "nrow", header.nrow);
  }

  // Source header read
  void read(XMLReader& xml, const string& path, PropSink_t& header)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
      /**************************************************************************/
    case 4:
      /**************************************************************************/
      break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "PropSink parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "wave_state", header.wave_state);
    if (header.wave_state != WAVE_TYPE_S_WAVE)
    {
      read(paramtop, "direction",  header.direction);
    }

    read(paramtop, "sink_type", header.sink_type);
    if (header.sink_type == SNK_TYPE_SHELL_SINK)
    {
      XMLReader shelltop(paramtop, "ShellSink");

      read(shelltop, "SinkSmearingParam", header.sinkSmearParam);
      read(shelltop, "laplace_power", header.laplace_power);
      read(shelltop, "link_smear_fact", header.link_smear_fact);
      read(shelltop, "link_smear_num", header.link_smear_num);
      read(shelltop, "disp_length", header.disp_length);
      read(shelltop, "disp_dir", header.disp_dir);
    }		
    read(paramtop, "nrow", header.nrow);
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
      read(paramtop, "nonRelSeqProp", param.nonRelSeqProp);
      read(paramtop, "seq_src", param.seq_src);
      break;

    default:
      QDPIO::cerr << "ChromaSeqProp parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "InvertParam", param.invParam);
    read(paramtop, "t_sink", param.t_sink);
    read(paramtop, "sink_mom", param.sink_mom);
    read(paramtop, "nrow", param.nrow);
  }


  //! SeqSource header reader
  void read(XMLReader& xml, const string& path, SeqSource_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
      /**************************************************************************/
    case 1:
      break;

    default:
      /**************************************************************************/
      QDPIO::cerr << "SeqSource parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "seq_src", param.seq_src);
    read(paramtop, "t_sink", param.t_sink);
    read(paramtop, "sink_mom", param.sink_mom);
    read(paramtop, "nrow", param.nrow);
  }


  // Forward propagator header read
  void read(XMLReader& xml, const string& path, ChromaProp_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    multi1d<int> boundary;

    switch (version) 
    {
      /**************************************************************************/
    case 5:
    {
      // In this modified version of v4, the fermion action specific stuff
      // goes into a <FermionAction> tag beneath <Param>
      param.nonRelProp = false;

      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);
      read(paramtop, "nrow", param.nrow);

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
      read(paramtop, "nonRelProp", param.nonRelProp); // new - is this prop non-relativistic
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);
      read(paramtop, "nrow", param.nrow);

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
      read(paramtop, "nrow", param.nrow);
      read(paramtop, "nonRelProp", param.nonRelProp);

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
  }


  // Forward propagator header read
  void read(XMLReader& xml, const string& path, ChromaMultiProp_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    multi1d<int> boundary;

    switch (version) 
    {
      /**************************************************************************/
    case 5:  // Backward compatibility with non Multi ChromaProp_t
      /**************************************************************************/
    {
      read(paramtop, "MultiMasses", param.MultiMasses);
    
      param.nonRelProp = false;
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);
      read(paramtop, "nrow", param.nrow);

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

      read(paramtop, "nonRelProp", param.nonRelProp); // new - is this prop non-relativistic
      read(paramtop, "InvertParam", param.invParam);
      read(paramtop, "boundary", boundary);
      read(paramtop, "nrow", param.nrow);

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
      read(paramtop, "nrow", param.nrow);
      read(paramtop, "nonRelProp", param.nonRelProp);

      if (paramtop.count("boundary") != 0)
      {
	QDPIO::cerr << "ChromaMultiProp: paranoia check - found a misplaced boundary" << endl; 
	QDP_abort(1);
      }
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




  //---------------------------------------------------------------------------
  // Source header writer
  void write(XMLWriter& xml, const string& path, const PropSource_t& header)
  {
    push(xml, path);

    int version = 5;
    write(xml, "version", version);
    write(xml, "wave_state", header.wave_state);
    if (header.wave_state != WAVE_TYPE_S_WAVE)
    {
      write(xml, "direction",  header.direction);
    }

    write(xml, "source_type", header.source_type);
    if (header.source_type == SRC_TYPE_SHELL_SOURCE)
    {
      push(xml, "ShellSource");
      write(xml, "SourceSmearingParam", header.sourceSmearParam);
      write(xml, "laplace_power", header.laplace_power);
      write(xml, "link_smear_fact", header.link_smear_fact);
      write(xml, "link_smear_num", header.link_smear_num);
      write(xml, "disp_length", header.disp_length);
      write(xml, "disp_dir", header.disp_dir);
      pop(xml);
    }

    write(xml, "j_decay",  header.j_decay);
    write(xml, "t_source",  header.t_source);

    write(xml, "nrow",  header.nrow);

    pop(xml);
  }


  // Source header writer
  void write(XMLWriter& xml, const string& path, const PropSink_t& header)
  {
    push(xml, path);

    int version = 4;
    write(xml, "version", version);
    write(xml, "wave_state", header.wave_state);
    if (header.wave_state != WAVE_TYPE_S_WAVE)
    {
      write(xml, "direction",  header.direction);
    }

    write(xml, "sink_type", header.sink_type);
    if (header.sink_type == SNK_TYPE_SHELL_SINK)
    {
      push(xml, "ShellSink");
      write(xml, "SinkSmearingParam", header.sinkSmearParam);
      write(xml, "laplace_power", header.laplace_power);
      write(xml, "link_smear_fact", header.link_smear_fact);
      write(xml, "link_smear_num", header.link_smear_num);
      write(xml, "disp_length", header.disp_length);
      write(xml, "disp_dir", header.disp_dir);
      pop(xml);
    }

    write(xml, "nrow",  header.nrow);

    pop(xml);
  }


  // Write propagator inversion parameters
  void write(XMLWriter& xml, const string& path, const ChromaProp_t& header)
  {
    push(xml, path);

    int version = 7;
    write(xml, "version", version);
    write(xml, "nonRelProp", header.nonRelProp); // new - is this prop non-relativistic
    xml << header.fermact;
    write(xml, "InvertParam", header.invParam);
    write(xml, "nrow", header.nrow);

    pop(xml);
  }

  // Write propagator inversion parameters
  void write(XMLWriter& xml, const string& path, const ChromaMultiProp_t& header)
  {
    push(xml, path);

    int version = 7;
    write(xml, "version", version);
    write(xml, "nonRelProp", header.nonRelProp); // new - is this prop non-relativistic
    write(xml, "MultiMasses", header.MultiMasses);
    xml << header.fermact;
    write(xml, "InvertParam", header.invParam);
    write(xml, "nrow", header.nrow);

    pop(xml);
  }


  //! SeqPropagator header writer
  void write(XMLWriter& xml, const string& path, const ChromaSeqProp_t& param)
  {
    push(xml, path);

    int version = 3;
    write(xml, "version", version);
    write(xml, "nonRelSeqProp", param.nonRelSeqProp);
    write(xml, "seq_src", param.seq_src);
    write(xml, "InvertParam", param.invParam);
    write(xml, "t_sink", param.t_sink);
    write(xml, "sink_mom", param.sink_mom);
    write(xml, "nrow", param.nrow);

    pop(xml);
  }


  //! SeqSource header writer
  void write(XMLWriter& xml, const string& path, const SeqSource_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "seq_src", param.seq_src);
    write(xml, "t_sink", param.t_sink);
    write(xml, "sink_mom", param.sink_mom);
    write(xml, "nrow", param.nrow);

    pop(xml);
  }


  //! SeqPropagator header writer
  void write(XMLWriter& xml, const string& path, const ForwardProp_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "PropSink", param.sink_header);
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

