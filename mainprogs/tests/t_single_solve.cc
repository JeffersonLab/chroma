#include "chroma.h"

// Program calculates just one spin/color position
// Input: like hadspec.ini.cml: Cfg, Source, Propagator

using namespace Chroma;

struct SingleSolve_input_t
{
  multi1d<int> nrow;
  GroupXML_t   cfg;
  QDP::Seed    rng_seed;
  GroupXML_t   source;
  GroupXML_t   fermact;
  GroupXML_t   invParam;
  int          color;
  int          spin;
};


bool linkageHack(void)
{
  bool foo = true;

  foo &= GaugeInitEnv::registerAll();
  foo &= QuarkSourceConstructionEnv::registerAll();
  foo &= WilsonTypeFermActsEnv::registerAll();

  return foo;
}


int findMeasurement(XMLReader& meas, const std::string& name)
{
  int n = meas.count("elem");
  for(int i = 1; i <= n; ++i)
  {
    std::ostringstream path;
    path << "elem[" << i << "]";
    XMLReader elem(meas, path.str());

    std::string elem_name;
    read(elem, "Name", elem_name);
    if (elem_name == name)
      return i;
  }

  throw std::string("No measurement with Name ") + name + " found in InlineMeasurements";
}


void read(XMLReader& xml, const std::string& path, SingleSolve_input_t& p)
{
  XMLReader top(xml, path);

  read(top, "Param/nrow", p.nrow);
  p.cfg = readXMLGroup(top, "Cfg", "cfg_type");

  if (top.count("RNG") > 0)
    read(top, "RNG", p.rng_seed);
  else
    p.rng_seed = 11;

  XMLReader meas(top, "Param/InlineMeasurements");

  {
    std::ostringstream path;
    path << "elem[" << findMeasurement(meas, "MAKE_SOURCE") << "]/Param";
    XMLReader src(meas, path.str());
    p.source = readXMLGroup(src, "Source", "SourceType");
  }

  {
    std::ostringstream path;
    path << "elem[" << findMeasurement(meas, "PROPAGATOR") << "]/Param";
    XMLReader prop(meas, path.str());
    p.fermact  = readXMLGroup(prop, "FermionAction", "FermAct");
    p.invParam = readXMLGroup(prop, "InvertParam", "invType");
  }

  p.color = 0;
  p.spin  = 0;
  if (top.count("SingleSolve") > 0)
  {
    XMLReader ss(top, "SingleSolve");
    if (ss.count("color") > 0)
      read(ss, "color", p.color);
    if (ss.count("spin") > 0)
      read(ss, "spin", p.spin);
  }

  if (p.color < 0 || p.color >= Nc || p.spin < 0 || p.spin >= Ns)
    throw std::string("SingleSolve color or spin out of range");
}


int main(int argc, char *argv[])
{
  Chroma::initialize(&argc, &argv);

  START_CODE();

  QDPIO::cout << "Linkage = " << linkageHack() << std::endl;

  SingleSolve_input_t input;
  XMLReader xml_in;
  try
  {
    xml_in.open(Chroma::getXMLInputFileName());
    read(xml_in, "/chroma", input);
  }
  catch(const std::string& e)
  {
    QDPIO::cerr << "T_SINGLE_SOLVE: Caught Exception reading XML: " << e << std::endl;
    QDP_abort(1);
  }

  Layout::setLattSize(input.nrow);
  Layout::create();

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "t_single_solve");
  write(xml_out, "Input", xml_in);
  proginfo(xml_out);

  QDP::RNG::setrn(input.rng_seed);
  write(xml_out, "RNG", input.rng_seed);

  multi1d<LatticeColorMatrix> u(Nd);
  try
  {
    XMLReader gauge_file_xml, gauge_xml;
    std::istringstream xml_c(input.cfg.xml);
    XMLReader cfgtop(xml_c);
    QDPIO::cout << "Gauge initialization: cfg_type = " << input.cfg.id << std::endl;

    Handle<GaugeInit> gaugeInit(TheGaugeInitFactory::Instance().createObject(input.cfg.id,
                                                                             cfgtop,
                                                                             input.cfg.path));
    (*gaugeInit)(gauge_file_xml, gauge_xml, u);
    write(xml_out, "Config_info", gauge_xml);
  }
  catch(const std::string& e)
  {
    QDPIO::cerr << "T_SINGLE_SOLVE: Caught Exception during gauge init: " << e << std::endl;
    QDP_abort(1);
  }

  MesPlq(xml_out, "Observables", u);

  LatticePropagator quark_source;
  try
  {
    std::istringstream xml_s(input.source.xml);
    XMLReader sourcetop(xml_s);
    QDPIO::cout << "Source = " << input.source.id << std::endl;

    Handle< QuarkSourceConstruction<LatticePropagator> >
      sourceConstruction(ThePropSourceConstructionFactory::Instance().createObject(input.source.id,
                                                                                   sourcetop,
                                                                                   input.source.path));
    quark_source = (*sourceConstruction)(u);
  }
  catch(const std::string& e)
  {
    QDPIO::cerr << "T_SINGLE_SOLVE: Caught Exception creating source: " << e << std::endl;
    QDP_abort(1);
  }

  LatticeFermion chi;
  PropToFerm(quark_source, chi, input.color, input.spin);

  Double norm_chi = sqrt(norm2(chi));
  QDPIO::cout << "Source color = " << input.color << "  spin = " << input.spin
              << "  |chi| = " << norm_chi << std::endl;
  if (toBool(norm_chi == Double(0)))
  {
    QDPIO::cerr << "T_SINGLE_SOLVE: source component is zero" << std::endl;
    QDP_abort(1);
  }
  chi *= Double(1) / norm_chi;

  LatticeFermion psi = zero;
  SystemSolverResults_t res;

  try
  {
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    std::istringstream xml_f(input.fermact.xml);
    XMLReader fermacttop(xml_f);
    QDPIO::cout << "FermAct = " << input.fermact.id << std::endl;

    Handle< FermionAction<T,P,Q> >
      S_f(TheFermionActionFactory::Instance().createObject(input.fermact.id,
                                                           fermacttop,
                                                           input.fermact.path));

    Handle< FermState<T,P,Q> > state(S_f->createState(u));

    StopWatch swatch;
    swatch.reset();
    swatch.start();
    Handle< SystemSolver<T> > solver(S_f->qprop(state, input.invParam));
    swatch.stop();
    QDPIO::cout << "Solver setup time = " << swatch.getTimeInSeconds() << " secs" << std::endl;
    double setup_time = swatch.getTimeInSeconds();

    swatch.reset();
    swatch.start();
    res = (*solver)(psi, chi);
    swatch.stop();
    QDPIO::cout << "Solve time = " << swatch.getTimeInSeconds() << " secs" << std::endl;

    push(xml_out, "Solve");
    write(xml_out, "color", input.color);
    write(xml_out, "spin", input.spin);
    write(xml_out, "n_count", res.n_count);
    write(xml_out, "resid", res.resid);
    write(xml_out, "setup_time", setup_time);
    write(xml_out, "solve_time", swatch.getTimeInSeconds());
    pop(xml_out);
  }
  catch(const std::string& e)
  {
    QDPIO::cerr << "T_SINGLE_SOLVE: Caught Exception during solve: " << e << std::endl;
    QDP_abort(1);
  }

  psi *= norm_chi;

  {
    SftMom phases(0, true, Nd-1);
    multi1d<Double> sol_corr = sumMulti(localNorm2(psi), phases.getSet());

    push(xml_out, "Solution_correlator");
    write(xml_out, "sol_corr", sol_corr);
    pop(xml_out);
  }

  QDPIO::cout << "n_count = " << res.n_count << "  resid = " << res.resid << std::endl;

  pop(xml_out);
  xml_out.close();

  END_CODE();

  Chroma::finalize();
  exit(0);
}
