#include "actions/ferm/invert/quda_solvers/quda_eig_params.h"

using namespace QDP;

namespace Chroma
{
  void read(XMLReader& xml_in, const std::string& path, QudaEigParam& p)
  {
    XMLReader paramtop(xml_in, path);
    readOptional<::QudaEigType>(paramtop, "./EigType", p.eig_type, QUDA_EIG_TR_LANCZOS);
    readQudaBool(paramtop, "./UsePolyAcc", p.use_poly_acc, true);
    readOptional(paramtop, "./PolyDeg", p.poly_deg,100);
    readOptional(paramtop, "./EigAmax", p.a_max, (double)4.0);
    readOptional(paramtop, "./EigAmin", p.a_min, (double)0.1);
    readQudaBool(paramtop, "./PreserveDeflation", p.preserve_deflation, false);
    readQudaBool(paramtop, "./PreserveEvals", p.preserve_evals, false);
    readQudaBool(paramtop, "./UseDagger", p.use_dagger, false);
    readQudaBool(paramtop, "./UseNormOp", p.use_norm_op, false);
    readQudaBool(paramtop, "./UsePC", p.use_pc, false);
    readQudaBool(paramtop, "./UseEigenQR", p.use_eigen_qr,false);
    readQudaBool(paramtop, "./ComputeSVD", p.compute_svd, false);
    readQudaBool(paramtop, "./ComputeGamma5", p.compute_gamma5, false);
    readQudaBool(paramtop, "./RequireConvergence", p.require_convergence, true);
    readOptional(paramtop, "./SpectrumType", p.spectrum, QUDA_SPECTRUM_SR_EIG); // Deflation uses SR mostly
    readOptional(paramtop, "./NEv", p.n_ev);
    readOptional(paramtop, "./NKr", p.n_kr );
    readOptional(paramtop, "./NLockedMax", p.nLockedMax);
    readOptional(paramtop, "./NConv", p.n_conv);
    readOptional(paramtop, "./NEvDeflate", p.n_ev_deflate, p.n_conv);
    readOptional(paramtop, "./Tol", p.tol, (double)(1.0e-6));
    readOptional(paramtop, "./QRTol", p.qr_tol,(double)(1.0e-11));
    readOptional(paramtop, "./CheckInterval", p.check_interval);
    readOptional(paramtop, "./MaxRestarts", p.max_restarts);
    readOptional(paramtop, "./BatchedRotate", p.batched_rotate);
    readOptional(paramtop, "./BlockSize", p.block_size,1);
    readOptional(paramtop, "./ComputeEvalsBatchSize", p.compute_evals_batch_size);
    readOptional(paramtop, "./MaxOrthoAttempts", p.max_ortho_attempts);
    readOptional(paramtop, "./OrthoBlockSize", p.ortho_block_size);
    QudaPrecisionType tmp_prec;
    readOptional(paramtop, "./SavePrecision", tmp_prec, DOUBLE);
    p.save_prec = theChromaToQudaPrecisionTypeMap::Instance()[tmp_prec];
    // Boilerplate things to set
    strcpy(p.vec_infile, "");
    strcpy(p.vec_outfile, "");
    p.io_parity_inflate = QUDA_BOOLEAN_FALSE;
    p.partfile = QUDA_BOOLEAN_FALSE; // ignored
    p.struct_size = sizeof(p);

  }

  void write(XMLWriter& xml, const std::string& path, const QudaEigParam& p)
  {
    push(xml, path);
    switch( p.eig_type) {
      case QUDA_EIG_TR_LANCZOS:
        write(xml, "EigType", "QUDA_EIG_TR_LANCZOS");
        break;
      case QUDA_EIG_BLK_TR_LANCZOS:
        write(xml, "EigType", "QUDA_EIG_BLK_TR_LANCZOS");
        break;
      case QUDA_EIG_IR_ARNOLDI:
        write(xml, "EigType", "QUDA_EIG_IR_ARNOLDI");
        break;
      case QUDA_EIG_BLK_IR_ARNOLDI:
        write(xml, "EigType", "QUDA_EIG_BLK_IR_ARNOLDI");
        break;
      default:
        break;
    };

    write(xml, "UsePolyAcc", p.use_poly_acc);
    write(xml, "PolyDeg", p.poly_deg);
    write(xml, "Amax", p.a_max);
    write(xml, "Amin", p.a_min);
    write(xml, "PreserveDeflation", p.preserve_deflation);
    write(xml, "PreserveEvals", p.preserve_evals);
    write(xml, "UseDagger", p.use_dagger);
    write(xml, "UseNormOp", p.use_norm_op);
    write(xml, "UsePC", p.use_pc);
    write(xml, "UseEigenQR", p.use_eigen_qr);
    write(xml, "ComputeSVD", p.compute_svd);
    write(xml, "ComputeGamma5", p.compute_gamma5);
    write(xml, "RequireConvergence", p.require_convergence);
    write(xml, "SpectrumType", p.spectrum);
    write(xml, "NEv", p.n_ev);
    write(xml, "NKr", p.n_kr );
    write(xml, "NLockedMax", p.nLockedMax);
    write(xml, "NConv", p.n_conv);
    write(xml, "NEvDeflate", p.n_ev_deflate);
    write(xml, "Tol", p.tol);
    write(xml, "QRTol", p.qr_tol);
    write(xml, "CheckInterval", p.check_interval);
    write(xml, "MaxRestarts", p.max_restarts);
    write(xml, "BatchedRotate", p.batched_rotate);
    write(xml, "BlockSize", p.block_size);
    write(xml, "ComputeEvalsBatchSize", p.compute_evals_batch_size);
    write(xml, "MaxOrthoAttempts", p.max_ortho_attempts);
    write(xml, "OrthoBlockSize", p.ortho_block_size);

    switch(p.save_prec) {
      case QUDA_QUARTER_PRECISION : 
        write(xml, "SavePrecision", "QUARTER");
        break;
      case QUDA_HALF_PRECISION : 
        write(xml, "SavePrecision", "HALF");
        break;
      case QUDA_SINGLE_PRECISION : 
        write(xml, "SavePrecision", "SINGLE");
        break;
      case QUDA_DOUBLE_PRECISION : 
        write(xml, "SavePrecision", "DOUBLE");
        break;

    };
    pop(xml);

  }

} // Namespace Chroma
