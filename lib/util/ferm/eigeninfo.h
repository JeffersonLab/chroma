#ifndef EIGENSTATE_H
#define EIGENSTATE_H


namespace Chroma
{
 
  class EigenInfo
    {
    public:
      EigenInfo() {}

      EigenInfo(multi1d<Real>& eval, Real& l, multi1d<LatticeFermion>& evec): evalues(eval), largest(l), evectors(evec) {
	if (eval.size() != evec.size() ) {
	  QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	  QDP_abort(1);
	}
      }

      ~EigenInfo() {}
      
      EigenInfo(const EigenInfo& e): evalues(e.evalues), largest(e.largest), evectors(e.evectors) {
	if (e.evalues.size() != e.evectors.size() ) {
	  QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	  QDP_abort(1);
	}
      }


      multi1d<Real>& getEvalues() {
	return evalues;
      }

      Real& getLargest() {
	return largest;
      }


      multi1d<LatticeFermion>& getEvectors() {
	return evectors;
      }

    private:
      multi1d<Real> evalues;
      multi1d<LatticeFermion> evectors;
      Real largest;
    };
}

















#endif
