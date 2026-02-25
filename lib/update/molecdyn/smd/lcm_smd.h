// -*- C++ -*-
/*! \file
 * \brief SMD trajectory
 *
 * Stochastic Molecular Dynamics trajectory
 */

#ifndef LAT_COL_MAT_SMD_NEW_H
#define LAT_COL_MAT_SMD_NEW_H

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/hmc/global_metropolis_accrej.h"
#include "actions/ferm/invert/mg_solver_exception.h"
#include "util/gauge/taproj.h"
#include "io/xmllog_io.h"
#include "handle.h"

namespace Chroma
{

  //! SMD trajectory
  /*! @ingroup molecdyn */
  class LatColMatSMDTrj
  {
  public:
    //! Constructor
    LatColMatSMDTrj(Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& _H_MC,
                    Handle< AbsMDIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& _MD_int,
                    const Real& _gamma,
                    const Real& _pf_gamma,
                    const InternalFieldsRefreshMode _pf_refresh_mode,
                    bool _accept_reject,
                    bool _measure_actions)
      : the_MD(_MD_int), the_H_MC(_H_MC), gamma(_gamma), pf_gamma(_pf_gamma),
        pf_refresh_mode(_pf_refresh_mode), internal_fields_initialized(false),
        accept_reject(_accept_reject),
        measure_actions(_measure_actions)
    {
      if (accept_reject && !measure_actions) {
        QDPIO::cerr << "SMD: AcceptReject requires MeasureActions" << std::endl;
        QDP_abort(1);
      }
    }

    //! Destructor
    ~LatColMatSMDTrj(void) {}

    //! Access the Hamiltonian used by this trajectory object
    AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& getMCHamiltonian(void)
    {
      return *the_H_MC;
    }

    //! Mark whether pseudofermion internal fields are already initialized
    void setInternalFieldsInitialized(bool value)
    {
      internal_fields_initialized = value;
    }

    //! Do the SMD trajectory
    void operator()(AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s,
                    const bool WarmUpP,
                    const bool CheckRevP)
    {
      START_CODE();
      StopWatch swatch;

      AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& MD = *the_MD;
      AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& H_MC = *the_H_MC;

      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      XMLWriter& xml_log = TheXMLLogWriter::Instance();

      push(xml_out, "SMDTrajectory");
      push(xml_log, "SMDTrajectory");

      write(xml_out, "WarmUpP", WarmUpP);
      write(xml_log, "WarmUpP", WarmUpP);

      // Refresh momenta using OU update
      swatch.reset();
      swatch.start();
      refreshP(s, MD.getTrajLength());
      swatch.stop();
      QDPIO::cout << "SMD_TIME: Momentum Refresh Time: " << swatch.getTimeInSeconds() << " \n";

      // Refresh pseudofermions
      bool internal_fields_pushed = false;
      try {
        swatch.reset();
        swatch.start();
        if (pf_refresh_mode == INTERNAL_FIELDS_REFRESH_OU && !internal_fields_initialized) {
          H_MC.refreshInternalFields(s);
        }
        else {
          H_MC.refreshInternalFields(s, MD.getTrajLength(), pf_gamma, pf_refresh_mode);
        }
        internal_fields_initialized = true;
        H_MC.pushInternalFields();
        internal_fields_pushed = true;
        swatch.stop();
        QDPIO::cout << "SMD_TIME: Pseudofermion Refresh Time: " << swatch.getTimeInSeconds() << " \n";
      }
      catch (MGSolverException e) {
        QDPIO::cout << "ERROR: Caught MG Solver exception in pseudofermion refresh" << std::endl;
        QDPIO::cout << "ERROR: Exception Was: " << e.whatStr() << std::endl;
        QDPIO::cout << "Aborting";
        QDP_abort(2);
      }

      Handle< AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > s_old(s.clone());

      bool acceptTraj = true;
      Double KE_old, PE_old;

      if (measure_actions) {
        try {
          push(xml_out, "H_old");
          push(xml_log, "H_old");
          swatch.reset();
          swatch.start();
          H_MC.mesE(*s_old, KE_old, PE_old);
          swatch.stop();
          QDPIO::cout << "SMD_TIME: Start Energy Time: " << swatch.getTimeInSeconds() << " \n";

          write(xml_out, "KE_old", KE_old);
          write(xml_log, "KE_old", KE_old);
          write(xml_out, "PE_old", PE_old);
          write(xml_log, "PE_old", PE_old);

          pop(xml_log);
          pop(xml_out);
        }
        catch (MGSolverException e) {
          QDPIO::cout << "ERROR: Caught MG Solver exception in Start Energy Calculation" << std::endl;
          QDPIO::cout << "ERROR: Exception Was: " << e.whatStr() << std::endl;
          QDPIO::cout << "Aborting";
          QDP_abort(2);
        }
      }

      try {
        swatch.start();
        MD.copyFields();
        MD(s, MD.getTrajLength());
        swatch.stop();
        QDPIO::cout << "SMD_TIME: Traj MD Time: " << swatch.getTimeInSeconds() << " \n";

        if (measure_actions) {
          swatch.reset();
          swatch.start();
          Double KE, PE;
          push(xml_out, "H_new");
          push(xml_log, "H_new");
          H_MC.mesE(s, KE, PE);
          write(xml_out, "KE_new", KE);
          write(xml_log, "KE_new", KE);
          write(xml_out, "PE_new", PE);
          write(xml_log, "PE_new", PE);
          pop(xml_log);
          pop(xml_out);
          swatch.stop();
          QDPIO::cout << "SMD_TIME: Finish Energy Time: " << swatch.getTimeInSeconds() << " \n";

          Double DeltaKE = KE - KE_old;
          Double DeltaPE = PE - PE_old;
          Double DeltaH  = DeltaKE + DeltaPE;
          Double AccProb = where(DeltaH < 0.0, Double(1), exp(-DeltaH));

          write(xml_out, "deltaKE", DeltaKE);
          write(xml_log, "deltaKE", DeltaKE);
          write(xml_out, "deltaPE", DeltaPE);
          write(xml_log, "deltaPE", DeltaPE);
          write(xml_out, "deltaH", DeltaH);
          write(xml_log, "deltaH", DeltaH);
          write(xml_out, "AccProb", AccProb);
          write(xml_log, "AccProb", AccProb);

          QDPIO::cout << "Delta H = " << DeltaH << std::endl;
          QDPIO::cout << "AccProb = " << AccProb << std::endl;

          if (!WarmUpP && accept_reject) {
            acceptTraj = globalMetropolisAcceptReject(DeltaH);
          }
        }
      }
      catch (MGSolverException e) {
        QDPIO::cout << "WARNING: Caught MG Solver Convergence Exception during MD or Final Energy Calculation!" << std::endl;
        QDPIO::cout << "WARNING: Exception was: " << e.whatStr() << std::endl;
        QDPIO::cout << "WARNING: Aborting" << std::endl;
        QDP_abort(2);
      }

      write(xml_out, "AcceptP", acceptTraj);
      write(xml_log, "AcceptP", acceptTraj);
      QDPIO::cout << "AcceptP = " << acceptTraj << std::endl;

      if (CheckRevP) {
        Handle< AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > s_rev(s.clone());

        try {
          swatch.reset();
          swatch.start();
          QDPIO::cout << "Reversing trajectory for reversability test" << std::endl;

          flipMomenta(*s_rev);
          MD(*s_rev, MD.getTrajLength());
          flipMomenta(*s_rev);

          Double KE_rev;
          Double PE_rev;
          H_MC.mesE(*s_rev, KE_rev, PE_rev);

          Double DeltaDeltaKE = KE_rev - KE_old;
          Double DeltaDeltaPE = PE_rev - PE_old;
          Double DeltaDeltaH = DeltaDeltaKE + DeltaDeltaPE;

          Double dq;
          Double dp;
          reverseCheckMetrics(dq, dp, *s_rev, *s_old);

          push(xml_log, "ReversibilityMetrics");
          write(xml_log, "DeltaDeltaH", fabs(DeltaDeltaH));
          write(xml_log, "DeltaDeltaKE", fabs(DeltaDeltaKE));
          write(xml_log, "DeltaDeltaPE", fabs(DeltaDeltaPE));
          write(xml_log, "DeltaQPerSite", dq);
          write(xml_log, "DeltaPPerSite", dp);
          pop(xml_log);

          QDPIO::cout << "Reversibility: DeltaDeltaH = " << fabs(DeltaDeltaH) << std::endl;
          QDPIO::cout << "Reversibility: DeltaQ      = " << dq << std::endl;
          QDPIO::cout << "Reversibility: DeltaP      = " << dp << std::endl;
          swatch.stop();
          QDPIO::cout << "SMD_TIME: Reverse Check Time: " << swatch.getTimeInSeconds() << " \n";
        }
        catch (MGSolverException e) {
          QDPIO::cout << "WARNING: Caught MG Solver Exception in Reverse Trajectory" << std::endl;
          QDPIO::cout << "WARNING: Exception was: " << e.whatStr() << std::endl;
          QDPIO::cout << "WARNING: Aborting" << std::endl;
          QDP_abort(2);
        }
      }

      if (!acceptTraj) {
        s.getQ() = s_old->getQ();
        s.getP() = s_old->getP();
        flipMomenta(s);
        if (internal_fields_pushed) {
          H_MC.popInternalFields();
        }
      }
      else {
        if (internal_fields_pushed) {
          H_MC.dropInternalFields();
        }
      }

      pop(xml_log);
      pop(xml_out);

      END_CODE();
    }

  private:
    Handle< AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > the_MD;
    Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > the_H_MC;

    Real gamma;
    Real pf_gamma;
    InternalFieldsRefreshMode pf_refresh_mode;
    bool internal_fields_initialized;
    bool accept_reject;
    bool measure_actions;

    void refreshP(AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s,
                  const Real& traj_length) const
    {
      const Real c1 = exp(-gamma * traj_length);
      const Real c2 = sqrt(Real(1) - c1 * c1);

      for (int mu = 0; mu < Nd; mu++) {
        LatticeColorMatrix eta;
        gaussian(eta);
        eta *= sqrt(Real(0.5));
        taproj(eta);

        s.getP()[mu] = c1 * s.getP()[mu] + c2 * eta;
        taproj(s.getP()[mu]);
      }
    }

    void flipMomenta(AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
    {
      multi1d<LatticeColorMatrix>& p = s.getP();
      for (int mu = 0; mu < Nd; mu++) {
        p[mu] *= Real(-1);
      }
    }

    void reverseCheckMetrics(Double& deltaQ, Double& deltaP,
                             const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s,
                             const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s_old) const
    {
      multi1d<LatticeColorMatrix> DelQ(Nd);
      multi1d<LatticeColorMatrix> DelP(Nd);
      for (int mu = 0; mu < Nd; mu++) {
        DelQ[mu] = s.getQ()[mu] - s_old.getQ()[mu];
        DelP[mu] = s.getP()[mu] - s_old.getP()[mu];
      }

      deltaQ = sqrt(norm2(DelQ[0]));
      deltaP = sqrt(norm2(DelP[0]));
      for (int mu = 1; mu < Nd; mu++) {
        deltaQ += sqrt(norm2(DelQ[mu]));
        deltaP += sqrt(norm2(DelP[mu]));
      }

      deltaQ /= (Double(Nd) * Double(Layout::vol()));
      deltaP /= (Double(Nd) * Double(Layout::vol()));
    }
  };

} // end namespace Chroma

#endif
