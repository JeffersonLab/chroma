#include "chromabase.h"
#include <exception>
#include <string>
#include <sstream>
namespace Chroma  {

  class MGSolverException : public std::exception 
  {
  public: 
    MGSolverException(const Real& Mass, const std::string& SubspaceId, int Iters, const Real& Rsd, const Real& RsdTarget) : 
      _Mass(Mass), _SubspaceId(SubspaceId), _Iters(Iters), _Rsd(Rsd), _RsdTarget(RsdTarget) {}

    
    virtual const char* what() const throw()
    {
      std::ostringstream message;
      message << "MultiGrid Exception. SubspaceId=" << _SubspaceId
	      << " Mass=" << _Mass
	      << " Iters=" << _Iters
	      << " RsdTarget=" << _RsdTarget 
	      << " Rsd=" << _Rsd << std::endl;

      std::string message_string = message.str();

      return message_string.c_str();

    }

    const Real&  getMass() const { return _Mass; }
    const std::string&  getSubspaceID() const { return _SubspaceId; }
    const int& getIters() const { return _Iters; }
    const Real& getRsd() const { return _Rsd; }
    const Real& getRsdTarget() const { return _RsdTarget; }

  private:
    Real _Mass; // quark mass
    std::string _SubspaceId; // the subspace
    int _Iters; // Number of iterations taken
    Real _Rsd;  // Actual Residuum 
    Real _RsdTarget; // Target Residuum

  };


}
