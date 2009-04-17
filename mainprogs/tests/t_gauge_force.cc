#include "chroma.h"
#include <iostream>

int main(int argc, char *argv[])
{
  // Initialise QDP
  Chroma::initialize(&argc, &argv);

  // Setup a small lattice
  const int nrow_arr[] = {4, 4, 4, 4};
  multi1d<int> nrow(Nd);
  nrow=nrow_arr;
  Layout::setLattSize(nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    Cfg_t foo; foo.cfg_type=CFG_TYPE_DISORDERED;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGIN";
    gaugeStartup(file_xml, config_xml, u, foo);
  }

  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_gauge_force");

  multi1d<LatticeColorMatrix> p(Nd);


  // Get Periodic Gauge Boundaries
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;

  Handle< GaugeBC<P,Q> > gbc(new PeriodicGaugeBC<P,Q>);
  Handle< CreateGaugeState<P,Q> > cgs(new CreateSimpleGaugeState<P,Q>(gbc));
  Real betaMC = Real(5.7);
  RectGaugeAct S_g_MC(cgs, betaMC);

  multi1d<LatticeColorMatrix> dsdu_1(Nd);
  multi1d<LatticeColorMatrix> dsdu_2(Nd);

  dsdu_1 = zero;
  dsdu_2 = zero;

  Handle< GaugeState<P,Q> > state1(S_g_MC.createState(u));
  S_g_MC.deriv(dsdu_1, state1);
  Double S1 = S_g_MC.S(state1);

  // Test gauge invariance
  LatticeColorMatrix g;
  rgauge(u, g);

  Handle< GaugeState<P,Q> > state2(S_g_MC.createState(u));
  S_g_MC.deriv(dsdu_2, state2);
  Double S2 = S_g_MC.S(state2);

  push(xml_out, "ForceDiff");
  write(xml_out, "S1", S1);
  write(xml_out, "S2", S2);
  write(xml_out, "S_diff", norm2(S1-S2));
  pop(xml_out);

  for(int mu=0; mu < Nd; mu++) 
  { 
//    taproj(dsdu_1[mu]);
//    taproj(dsdu_2[mu]);

//    push(xml_out, "dsdu");
//    write(xml_out, "dsdu_1", dsdu_1[mu]);
//    write(xml_out, "dsdu_2", dsdu_2[mu]);
//    pop(xml_out);
  
    Double sum_diff=norm2(dsdu_1[mu] - adj(g)*dsdu_2[mu]*g);
    QDPIO::cout << "Mu = " << mu << " Sum Diff=" << sum_diff << endl;

    push(xml_out, "ForceDiff");
    write(xml_out, "mu", mu);
    write(xml_out, "dsdu1_norm", norm2(dsdu_1[mu]));
    write(xml_out, "dsdu2_norm", norm2(dsdu_2[mu]));
    write(xml_out, "dsdu_diff", sum_diff);
    pop(xml_out);
  }

  pop(xml_out);
  xml_out.close();

    // Finish
  Chroma::finalize();

  exit(0);
}

