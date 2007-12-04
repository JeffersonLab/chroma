// $Id: simple_baryon_2pt_w.cc,v 1.3 2007-12-04 02:59:05 kostas Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#include "meas/hadron/simple_baryon_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "meas/hadron/barspinmat_w.h"
#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/baryon_spinmat_funcmap_w.h"
#include "util/ft/sftmom.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SimpleBaryonSeqSourceEnv::Params& param)
  {
    SimpleBaryonSeqSourceEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const SimpleBaryonSeqSourceEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  namespace SimpleBaryonSeqSourceEnv
  { 
    //! The T and Spin struct
    /*! \ingroup hadron */
    struct SpinMatTsp_t
    {
      GroupXML_t  T_xml;          /*!< Holds source xml params*/
      SpinMatrix  T;

      GroupXML_t  sp_xml;         /*!< Holds source xml params*/
      SpinMatrix  sp;
    };
  }

  namespace GeneralBaryonSeqSourceEnv
  { 
    //! The T and Spin struct
    /*! \ingroup hadron */
    struct SpinMatTsp_t
    {
      GroupXML_t  T_xml;          /*!< Holds source xml params*/
      SpinMatrix  T;
      
      GroupXML_t  SRC_sp_xml;    /*!< Holds source xml params*/
      SpinMatrix  SRC_sp;

      GroupXML_t  SNK_sp_xml;    /*!< Holds source xml params*/
      SpinMatrix  SNK_sp;
    };
  }


  //! Read a T and sp struct
  /*! \ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleBaryonSeqSourceEnv::SpinMatTsp_t& param)
  {
    XMLReader paramtop(xml, path);

    param.T_xml   = readXMLGroup(paramtop, "Projector", "name");
    param.sp_xml  = readXMLGroup(paramtop, "DiquarkSpin", "name");

    param.T = BaryonSpinMatrixEnv::TheBarSpinMatFuncMap::Instance().callFunction(
      param.T_xml.id, paramtop, param.T_xml.path);
    param.sp = BaryonSpinMatrixEnv::TheBarSpinMatFuncMap::Instance().callFunction(
      param.sp_xml.id, paramtop, param.sp_xml.path);
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const SimpleBaryonSeqSourceEnv::SpinMatTsp_t& param)
  {
    push(xml, path);

    xml << param.T_xml.xml;
    xml << param.sp_xml.xml;

    pop(xml);
  }


  //! Read a T and sp struct
  /*! \ingroup hadron */
  void read(XMLReader& xml, const string& path, GeneralBaryonSeqSourceEnv::SpinMatTsp_t& param)
  {
    XMLReader paramtop(xml, path);

    param.T_xml   = readXMLGroup(paramtop, "Projector", "name");
    param.SRC_sp_xml  = readXMLGroup(paramtop, "SrcDiquarkSpin", "name");
    param.SRC_sp_xml  = readXMLGroup(paramtop, "SnkDiquarkSpin", "name");

    param.T = BaryonSpinMatrixEnv::TheBarSpinMatFuncMap::Instance().callFunction(param.T_xml.id, paramtop, param.T_xml.path);
    param.SRC_sp = BaryonSpinMatrixEnv::TheBarSpinMatFuncMap::Instance().callFunction(param.SRC_sp_xml.id, paramtop, param.SRC_sp_xml.path);
    param.SNK_sp = BaryonSpinMatrixEnv::TheBarSpinMatFuncMap::Instance().callFunction(param.SNK_sp_xml.id, paramtop, param.SNK_sp_xml.path);
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const SimpleBaryonSeqSourceEnv::SpinMatTsp_t& param)
  {
    push(xml, path);

    xml << param.T_xml.xml;
    xml << param.SRC_sp_xml.xml;
    xml << param.SNK_sp_xml.xml;

    pop(xml);
  }



  // Anonymous namespace
  /*! \ingroup hadron */
  namespace
  {
    //! Check only 1 prop passed
    /*! \ingroup hadron */
    void check1Args(const char* name, const multi1d<LatticePropagator>& quark_propagators)
    {
      if (quark_propagators.size() != 1)
      {
//      QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
//      QDP_abort(1);

	ostringstream s;
	s << name << ": expecting 1 prop, instead passed = " << quark_propagators.size() << endl;
	throw s.str();
      }
    }

    //! Check only 2 props passed
    /*! \ingroup hadron */
    void check2Args(const char* name, const multi1d<LatticePropagator>& quark_propagators)
    {
      if (quark_propagators.size() != 2)
      {
//      QDPIO::cerr << __func__ << ": expect only 2 prop" << endl;
//      QDP_abort(1);

	ostringstream s;
	s << name << ": expecting 2 props, instead passed = " << quark_propagators.size() << endl;
	throw s.str();
      }
    }
  }


  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace SimpleBaryonSeqSourceEnv
  { 

    //! Initialize
    Params::Params()
    {
      j_decay = -1;
      t_sink  = -1;
      sink_mom.resize(Nd-1);
      sink_mom = 0;
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "time_rev", param.time_rev);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_sink", t_sink);
      read(paramtop, "sink_mom", sink_mom);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "time_rev", param.time_rev);
      write(xml, "j_decay", j_decay);
      write(xml, "t_sink", t_sink);
      write(xml, "sink_mom", sink_mom);
      pop(xml);
    }


    // Nucleon-Nucleon 2pt correlator with general projector and Cg5
    multi1d<Hadron2PtContraction_t> operator()(const multi1d<LatticeColorMatrix>& u);
    {
      START_CODE();

      multi1d<ForwardProp_t>& forward_headers;
      multi1d<LatticePropagator>& quark_propagators;

      check2Args("BarNuclTCg5", quark_propagators);

      LatticePropagator src_prop_tmp;
      LatticePropagator q1_tmp;
      LatticePropagator q2_tmp;
      LatticePropagator di_quark;
      LatticeColorMatrix col_mat;
  
      /* "\bar u O u" insertion in NR proton, ie. 
       * "(u Cg5 d) u" */
      /* Some generic T */

      // Use precomputed Cg5
      q1_tmp = quark_propagators[0] * Cg5;
      q2_tmp = Cg5 * quark_propagators[1];
      di_quark = quarkContract24(q1_tmp, q2_tmp);

      // First term
      src_prop_tmp = T * di_quark;

      // Now the second term
      src_prop_tmp += traceSpin(di_quark) * T;

      // The third term...
      q1_tmp = q2_tmp * Cg5;
      q2_tmp = quark_propagators[0] * T;

      src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

      END_CODE();

      return projectBaryon(src_prop_tmp, forward_headers);
    }

    
    // Compute the 2-pt at the sink
    Complex
    BarNuclUTCg5::twoPtSink(const multi1d<LatticeColorMatrix>& u,
			    const multi1d<ForwardProp_t>& forward_headers,
			    const multi1d<LatticePropagator>& quark_propagators,
			    int gamma_insertion)
    {
      check2Args("BarNuclUTCg5", quark_propagators);
      setTSrce(forward_headers);
      setBC(forward_headers);

      // Constructor the 2pt
      LatticeComplex b_prop = Baryon2PtContractions::sigma2pt(quark_propagators[0], 
							      quark_propagators[1], 
							      T, Cg5);

      // Constructor the 2pt
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(b_prop);

      // U-quark rescaling factor is 1 for this type of seqsource
      return Real(2) * hsum[0][getTSink()];
    }


    //! Nucleon-Nucleon D piece with general projector and Cg5
    LatticePropagator
    BarNuclDTCg5::operator()(const multi1d<LatticeColorMatrix>& u,
			     const multi1d<ForwardProp_t>& forward_headers,
			     const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      check1Args("BarNuclDTCg5", quark_propagators);

      LatticePropagator src_prop_tmp;
      LatticePropagator q1_tmp;
      LatticePropagator q2_tmp;

      /* "\bar d O d" insertion in NR proton, ie. 
       * "(u Cg5 d) u" */
      /* Some generic T */

      // First term
      q2_tmp = quark_propagators[0] * Cg5;
      q1_tmp = T * q2_tmp;

      q2_tmp = Cg5 * quark_propagators[0];
      src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

      // Second term
      q1_tmp = q2_tmp * Cg5;
      q2_tmp = quark_propagators[0] * T;

      src_prop_tmp -= transposeSpin(quarkContract12(q2_tmp, q1_tmp));

      END_CODE();

      return projectBaryon(src_prop_tmp,
			   forward_headers);
    }

    // Compute the 2-pt at the sink
    Complex
    BarNuclDTCg5::twoPtSink(const multi1d<LatticeColorMatrix>& u,
			    const multi1d<ForwardProp_t>& forward_headers,
			    const multi1d<LatticePropagator>& quark_propagators,
			    int gamma_insertion)
    {
      check2Args("BarNuclDTCg5", quark_propagators);
      setTSrce(forward_headers);
      setBC(forward_headers);

      // Constructor the 2pt
      LatticeComplex b_prop = Baryon2PtContractions::sigma2pt(quark_propagators[0], 
							      quark_propagators[1], 
							      T, Cg5);

      // Extract the sink at the appropriate momenta
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(b_prop);

      // D-quark rescaling factor is 1 for this type of seqsource
      return Real(1) * hsum[0][getTSink()];
    }



    //! Delta+ - Delta+ U piece with general projector and spin matrix
    LatticePropagator
    BarDeltaUTsp::operator()(const multi1d<LatticeColorMatrix>& u,
			     const multi1d<ForwardProp_t>& forward_headers,
			     const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      check2Args("BarDeltaUTsp", quark_propagators);

      LatticePropagator src_prop_tmp;
      LatticePropagator q1_tmp;
      LatticePropagator q2_tmp;
      LatticePropagator di_quark;
      LatticeColorMatrix col_mat;
  
      /* "\bar u O u" insertion in Delta^+,
	 ie. "2*(u sp d) u + (u sp u) d" */
      // generic T

      q1_tmp = quark_propagators[0] * sp;
      q2_tmp = sp * quark_propagators[1];
      di_quark = quarkContract24(q1_tmp, q2_tmp);
      src_prop_tmp = T * di_quark;

      src_prop_tmp += traceSpin(di_quark) * T;

      q1_tmp = q2_tmp * sp;
      q2_tmp = quark_propagators[0] * T;
      src_prop_tmp += quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

      q1_tmp = sp * quark_propagators[0];
      q2_tmp = q1_tmp * sp;
      q1_tmp = T * quark_propagators[1];
      src_prop_tmp += transposeSpin(quarkContract12(q1_tmp, q2_tmp)) + quarkContract24(q1_tmp, q2_tmp);

      q2_tmp = q1_tmp * sp;
      q1_tmp = sp * quark_propagators[0];
      src_prop_tmp += quarkContract14(q2_tmp, q1_tmp);

      q1_tmp = sp * quark_propagators[1];
      q2_tmp = q1_tmp * T;
      q1_tmp = quark_propagators[0] * sp;
      src_prop_tmp += quarkContract14(q1_tmp, q2_tmp) + quarkContract13(q1_tmp, q2_tmp);
      src_prop_tmp *= 2;

      END_CODE();

      return projectBaryon(src_prop_tmp,
			   forward_headers);
    }


    // Compute the 2-pt at the sink
    Complex
    BarDeltaUTsp::twoPtSink(const multi1d<LatticeColorMatrix>& u,
			    const multi1d<ForwardProp_t>& forward_headers,
			    const multi1d<LatticePropagator>& quark_propagators,
			    int gamma_insertion)
    {
      check2Args("BarDeltaUTsp", quark_propagators);
      setTSrce(forward_headers);
      setBC(forward_headers);

      // Constructor the 2pt
      LatticeComplex b_prop = Baryon2PtContractions::sigmast2pt(quark_propagators[0], 
								quark_propagators[1], 
								T, sp);

      // Extract the sink at the appropriate momenta
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(b_prop);

      // U-quark rescaling factor is 1 for this type of seqsource
      return Real(2) * hsum[0][getTSink()];
    }



    //! Delta+ - Delta+ D piece with general projector and spin matrix
    LatticePropagator
    BarDeltaDTsp::operator()(const multi1d<LatticeColorMatrix>& u,
			     const multi1d<ForwardProp_t>& forward_headers,
			     const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      check1Args("BarDeltaDTsp", quark_propagators);

      LatticePropagator src_prop_tmp;
      LatticePropagator q1_tmp;
      LatticePropagator q2_tmp;
      LatticePropagator di_quark;

      /* "\bar d O d" insertion in Delta^+,
	 ie. "2*(u sp d) u + (u sp u) d" */
      // generic T

      q2_tmp = quark_propagators[0] * sp;
      q1_tmp = T * q2_tmp;
      q2_tmp = sp * quark_propagators[0];
      src_prop_tmp = quarkContract14(q1_tmp, q2_tmp);

      q1_tmp = q2_tmp * sp;
      q2_tmp = quark_propagators[0] * T;
      src_prop_tmp += transposeSpin(quarkContract12(q2_tmp, q1_tmp));

      q1_tmp = quark_propagators[0] * sp;
      q2_tmp = sp * quark_propagators[0];
      src_prop_tmp += T * quarkContract24(q1_tmp, q2_tmp);

      di_quark = quarkContract13(q1_tmp, q2_tmp);
      src_prop_tmp += di_quark * T;
      src_prop_tmp *= 2;

      src_prop_tmp += traceSpin(di_quark) * T;

      END_CODE();

      return projectBaryon(src_prop_tmp,
			   forward_headers);
    }

    // Compute the 2-pt at the sink
    Complex
    BarDeltaDTsp::twoPtSink(const multi1d<LatticeColorMatrix>& u,
			    const multi1d<ForwardProp_t>& forward_headers,
			    const multi1d<LatticePropagator>& quark_propagators,
			    int gamma_insertion)
    {
      check2Args("BarDeltaDTsp", quark_propagators);
      setTSrce(forward_headers);
      setBC(forward_headers);

      // Constructor the 2pt
      LatticeComplex b_prop = Baryon2PtContractions::sigmast2pt(quark_propagators[0], 
								quark_propagators[1], 
								T, sp);

      // Extract the sink at the appropriate momenta
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(b_prop);

      // D-quark rescaling factor is 1 for this type of seqsource
      return Real(1) * hsum[0][getTSink()];
    }



    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! "\bar u O u" insertion in nucleon-nucleon */
      /*!
       * \ingroup hadron
       *
       * This is a generic version
       */
      HadronSeqSource<LatticePropagator>* barNuclNuclU(XMLReader& xml_in,
						       const std::string& path)
      {
	// Determine the spin matrices
	SpinMatTsp_t spin;
	read(xml_in, path, spin);
    
	return new BarNuclUTCg5(Params(xml_in, path), spin.T, spin.sp);
      }


      //! "\bar d O d" insertion in nucleon-nucleon */
      /*!
       * \ingroup hadron
       *
       * This is a generic version
       */
      HadronSeqSource<LatticePropagator>* barNuclNuclD(XMLReader& xml_in,
						       const std::string& path)
      {
	// Determine the spin matrices
	SpinMatTsp_t spin;
	read(xml_in, path, spin);
    
	return new BarNuclDTCg5(Params(xml_in, path), spin.T, spin.sp);
      }


      //! "\bar u O u" insertion in proton, ie.  "(u C gamma_5 d) u" */
      /*!
       * \ingroup hadron
       *
       * C gamma_5 = Gamma(5) = - (C gamma_5)^T
       * T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 
       */
      HadronSeqSource<LatticePropagator>* barNuclUUnpol(XMLReader& xml_in,
							const std::string& path)
      {
	return new BarNuclUTCg5(Params(xml_in, path), BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5());
      }


      //! "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u"
      /*!
       * \ingroup hadron
       *
       * C gamma_5 = Gamma(5) = - (C gamma_5)^T
       * T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 
       */
      HadronSeqSource<LatticePropagator>* barNuclDUnpol(XMLReader& xml_in,
							const std::string& path)
      {
	return new BarNuclDTCg5(Params(xml_in, path), BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5());
      }


      //! "\bar u O u" insertion in proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       *
       * C g_5 = C gamma_5
       * T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2
       */
      HadronSeqSource<LatticePropagator>* barNuclUPol(XMLReader& xml_in,
						      const std::string& path)
      {
	return new BarNuclUTCg5(Params(xml_in, path), BaryonSpinMats::Tpol(), BaryonSpinMats::Cg5());
      }


      //! "\bar u O u" insertion in proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       *
       * C g_5 = C gamma_5
       * T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 
       */
      HadronSeqSource<LatticePropagator>* barNuclDPol(XMLReader& xml_in,
						      const std::string& path)
      {
	return new BarNuclDTCg5(Params(xml_in, path), BaryonSpinMats::Tpol(), BaryonSpinMats::Cg5());
      }

      //! "\bar u O u" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       *
       * C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
       * T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 
       */
      HadronSeqSource<LatticePropagator>* barNuclUUnpolNR(XMLReader& xml_in,
							  const std::string& path)
      {
	return new BarNuclUTCg5(Params(xml_in, path), BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5NR());
      }


      //! "\bar d O d" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       *
       * C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
       * T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 
       */
      HadronSeqSource<LatticePropagator>* barNuclDUnpolNR(XMLReader& xml_in,
							  const std::string& path)
      {
	return new BarNuclDTCg5(Params(xml_in, path), BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5NR());
      }


      //! "\bar u O u" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       *
       * C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
       * T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 
       */
      HadronSeqSource<LatticePropagator>* barNuclUPolNR(XMLReader& xml_in,
							const std::string& path)
      {
	return new BarNuclUTCg5(Params(xml_in, path), BaryonSpinMats::Tpol(), BaryonSpinMats::Cg5NR());
      }


      //! "\bar d O d" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       *
       * C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
       * T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 
       */
      HadronSeqSource<LatticePropagator>* barNuclDPolNR(XMLReader& xml_in,
							const std::string& path)
      {
	return new BarNuclDTCg5(Params(xml_in, path), BaryonSpinMats::Tpol(), BaryonSpinMats::Cg5NR());
      }


      //! \bar u O u" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" 
      /*!
       * \ingroup hadron
       * 
       * \f$C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )\f$
       * 
       * \f$T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       *      = (1 + Gamma(8) - i G(3) - i G(11)) / 2\f$
       */
      HadronSeqSource<LatticePropagator>* barNuclUMixedNR(XMLReader& xml_in,
							  const std::string& path)
      {
	return new BarNuclUTCg5(Params(xml_in, path), BaryonSpinMats::Tmixed(), BaryonSpinMats::Cg5NR());
      }


      //! "\bar d O d" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       * 
       * C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
       *
       * T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       *   = (1 + Gamma(8) - i G(3) - i G(11)) / 2
       */
      HadronSeqSource<LatticePropagator>* barNuclDMixedNR(XMLReader& xml_in,
							  const std::string& path)
      {
	return new BarNuclDTCg5(Params(xml_in, path), BaryonSpinMats::Tmixed(), BaryonSpinMats::Cg5NR());
      }

      //! \bar u O u" insertion in NR proton
      /*!
       * \ingroup hadron
       * 
       * "\bar u O u" insertion in NR proton, ie. 
       * "(u C gamma_5 (1/2)(1 - gamma_4)  d) u" 
       * 
       * \f$C g_5 NR = (1/2)*C gamma_5 * ( 1 - g_4 )\f$
       * 
       * $T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       *   = (1 + Gamma(8) - i G(3) - i G(11)) / 2$
       */
      HadronSeqSource<LatticePropagator>* barNuclUMixedNRnegPar(XMLReader& xml_in,
								const std::string& path)
      {
	return new BarNuclUTCg5(Params(xml_in, path), 
				BaryonSpinMats::TmixedNegPar(), BaryonSpinMats::Cg5NRnegPar());
      }

      //! "\bar d O d" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 - gamma_4)  d) u"
      /*!
       * \ingroup hadron
       * 
       * C g_5 NR = (1/2)*C gamma_5 * ( 1 - g_4 )
       *
       * T = (1 + \Sigma_3)*(1 - gamma_4) / 2 
       *   = (1 - Gamma(8) + i G(3) - i G(11)) / 2
       */
      HadronSeqSource<LatticePropagator>* barNuclDMixedNRnegPar(XMLReader& xml_in,
								const std::string& path)
      {
	return new BarNuclDTCg5(Params(xml_in, path), 
				BaryonSpinMats::TmixedNegPar(), BaryonSpinMats::Cg5NRnegPar());
      }


      /** The octet baryon sequantial sources **/
      //! \bar d O d" insertion in NR proton
      /*!
       * \ingroup hadron
       * 
       * "\bar d O d" insertion in NR proton, ie. 
       * "(s C gamma_5 (1/2)(1 + gamma_4)  d) s" 
       * 
       * $C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )$
       * 
       * $T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       *   = (1 + Gamma(8) - i G(3) - i G(11)) / 2$
   
       * The d-quark insertion for a Xi baryon: primarily for the
       * Xi- to Xi0 transision 

       * The Xi is just like the proton with up quark replaced with the strange
       * the single quark propagator passed in is just the strange quark propagator
       */
      HadronSeqSource<LatticePropagator>* barXiDMixedNR(XMLReader& xml_in,
							const std::string& path)
      {
	return new BarNuclDTCg5(Params(xml_in, path), BaryonSpinMats::Tmixed(), BaryonSpinMats::Cg5NR());
      }
  


      //! "\bar u O u" insertion in delta-delta */
      /*!
       * \ingroup hadron
       *
       * This is a generic version
       */
      HadronSeqSource<LatticePropagator>* barDeltaDeltaU(XMLReader& xml_in,
							 const std::string& path)
      {
	// Determine the spin matrices
	SpinMatTsp_t spin;
	read(xml_in, path, spin);
    
	return new BarDeltaUTsp(Params(xml_in, path), spin.T, spin.sp);
      }


      //! "\bar d O d" insertion in delta-delta */
      /*!
       * \ingroup hadron
       *
       * This is a generic version
       */
      HadronSeqSource<LatticePropagator>* barDeltaDeltaD(XMLReader& xml_in,
							 const std::string& path)
      {
	// Determine the spin matrices
	SpinMatTsp_t spin;
	read(xml_in, path, spin);
    
	return new BarDeltaDTsp(Params(xml_in, path), spin.T, spin.sp);
      }


      //! "\bar u O u" insertion in Delta^+ "2*(u sp d) u + (u sp u) d" */
      /*!
       * \ingroup hadron
       * 
       * C gamma_- = sp = (C gamma_-)^T
       * T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2
       *
       * Agghh, we have a goofy factor of 2 normalization factor here. The
       * ancient szin way didn't care about norms, so it happily made it
       * 2 times too big. There is a missing 0.5 in T_unpol guy.
       * Since nobody has used this code before, we are switching to a more
       * sane convention and breaking agreement with the old szin code.
       */
      HadronSeqSource<LatticePropagator>* barDeltaUUnpol(XMLReader& xml_in,
							 const std::string& path)
      {
	return new BarDeltaUTsp(Params(xml_in, path), BaryonSpinMats::Tunpol(), BaryonSpinMats::Cgm());
      }


      //! "\bar d O d" insertion in Delta^+ "2*(u sp d) u + (u sp u) d" */
      /*!
       * \ingroup hadron
       * 
       * C gamma_- = sp = (C gamma_-)^T 
       * T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2
       *
       *
       * Agghh, we have a goofy factor of 2 normalization factor here. The
       * ancient szin way didn't care about norms, so it happily made it
       * 2 times too big. There is a missing 0.5 in T_unpol guy.
       * Since nobody has used this code before, we are switching to a more
       * sane convention and breaking agreement with the old szin code.
       */
      HadronSeqSource<LatticePropagator>* barDeltaDUnpol(XMLReader& xml_in,
							 const std::string& path)
      {
	return new BarDeltaDTsp(Params(xml_in, path), BaryonSpinMats::Tunpol(), BaryonSpinMats::Cgm());
      }


      //! Local registration flag
      bool registered = false;

    }  // end anonymous namespace


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register needed stuff
	success &= BaryonSpinMatrixEnv::registerAll();

	//! Register all the factories
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL-NUCL_U"), 
										      barNuclNuclU);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL-NUCL_D"), 
										      barNuclNuclD);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_U_UNPOL"), 
										      barNuclUUnpol);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_D_UNPOL"), 
										      barNuclDUnpol);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_U_POL"),
										      barNuclUPol);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_D_POL"),
										      barNuclDPol);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_U_UNPOL_NONREL"),
										      barNuclUUnpolNR);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_D_UNPOL_NONREL"),
										      barNuclDUnpolNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_U_POL_NONREL"),
										      barNuclUPolNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_D_POL_NONREL"),
										      barNuclDPolNR);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_U_MIXED_NONREL"),
										      barNuclUMixedNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_D_MIXED_NONREL"),  
										      barNuclDMixedNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_U_MIXED_NONREL_NEGPAR"),
										      barNuclUMixedNRnegPar);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("NUCL_D_MIXED_NONREL_NEGPAR"),
										      barNuclDMixedNRnegPar);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("XI_D_MIXED_NONREL"),   
										      barXiDMixedNR);

      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("DELTA-DELTA_U"), 
										      barDeltaDeltaU);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("DELTA-DELTA_D"), 
										      barDeltaDeltaD);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("DELTA_U_UNPOL"),
										      barDeltaUUnpol);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("DELTA_D_UNPOL"),
										      barDeltaDUnpol);

	registered = true;
      }
      return success;
    }

  } // namespace BaryonSeqSourceCallMapEnv


}  // end namespace Chroma
