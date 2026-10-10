/*! \file
 * \brief Inline task to execute commands
 */

#include "chroma_config.h"

#ifdef BUILD_SB
// Activate the MPI support in Superbblas
#  define SUPERBBLAS_USE_MPI

// Activate redstar-datalib support for superbblas and with gpus
#  define USE_SUPERBBLAS
#  define USE_SUPERBBLAS_WITH_GPU_SUPPORT
#endif

#ifdef BUILD_REDSTAR_DATALIB
#  include "algs/superb_tensor.h"
#endif

#include "meas/inline/io/inline_cmd.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "qdp.h"
#include <cstdlib>

#include "util/ferm/superb_contractions.h"

namespace Chroma 
{ 
  namespace InlineExecuteCmdEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "CMD";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }

    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "cmd", input.cmd);
      input.only_on_master = false;
      if (inputtop.count("on_on_master") == 1)
      {
	read(inputtop, "only_on_master", input.only_on_master);
      }
      input.max_attempts = 1;
      if (inputtop.count("max_attempts") == 1)
      {
	read(inputtop, "max_attempts", input.max_attempts);
	if (input.max_attempts < 1)
	  throw std::runtime_error("invalid value of max_attempts");
      }
      input.constrain_to_gpu = false;
      if (inputtop.count("constrain_to_gpu") == 1)
      {
	read(inputtop, "constrain_to_gpu", input.constrain_to_gpu);
      }
    }

    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "cmd", input.cmd);
      write(xml, "only_on_master", input.only_on_master);
      write(xml, "max_attempts", input.max_attempts);
      write(xml, "constrain_to_gpu", input.constrain_to_gpu);

      pop(xml);
    }

    void read(XMLReader& xml, const std::string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    void write(XMLWriter& xml, const std::string& path, const Params& input)
    {
      push(xml, path);

      write(xml, "Param", input.param);

      pop(xml);
    }

    int execute(const std::string& cmd0, bool constrain_to_gpu)
    {
      std::string cmd;
      if (constrain_to_gpu)
      {
#if defined(SUPERBBLAS_USE_GPU)
	const int gpu_device = SB::detail::getGpuContext()->device;
	cmd =
#  ifdef SUPERBBLAS_USE_CUDA
	  std::string("CUDA_VISIBLE_DEVICES=") +
#  else
	  std::string("ROCR_VISIBLE_DEVICES=") +
#  endif
	  std::to_string(gpu_device) + std::string(" ") + cmd0;
#else
	throw std::runtime_error(
	  "chroma compiled without gpu support but requested constrain_to_gpu=true");
#endif // SUPERBBLAS_USE_GPU
      }
      else
      {
	cmd = cmd0;
      }
      std::cout << "on proc " << Layout::nodeNumber() << " executing " << cmd << std::endl;
      int ret = std::system(cmd.c_str());
#if defined(_WIN32)
      const auto exit_code = ret;
#else
      const auto exit_code = WEXITSTATUS(ret);
#endif
      std::cout << "on proc " << Layout::nodeNumber() << " result " << exit_code << std::endl;
      return exit_code;
    }

    // Param stuff
    Params::Params()
    {
      frequency = 0;
    }

    Params::Params(XMLReader& xml_in, const std::string& path)
    {
      try
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for source construction
	read(paramtop, "Param", param);
      } catch (const std::string& e)
      {
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }

    void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
    {
      START_CODE();

      QDPIO::cout << InlineExecuteCmdEnv::name << ": executing command";
      if (params.param.only_on_master)
	QDPIO::cout << " on master process";
      else
	QDPIO::cout << " one line on each process";
      QDPIO::cout << ":" << std::endl;

      if (params.param.only_on_master)
      {
	if (Layout::nodeNumber() == 0)
	{
	  for (const auto& cmd : params.param.cmd)
	  {
	    for (int attempt = 0; attempt < params.param.max_attempts; ++attempt)
	    {
	      if (execute(cmd, params.param.constrain_to_gpu) == 0)
	      {
		break;
	      }
	      else if (attempt == params.param.max_attempts - 1)
	      {
		throw std::runtime_error("some tasks failed!");
	      }
	    }
	  }
	}
      }
      else
      {
	std::vector<int> remaining_cmds;
	for (int i = 0; i < params.param.cmd.size(); ++i)
	  remaining_cmds.push_back(i);
	std::vector<int> procs;
	for (int attempt = 0; attempt < params.param.max_attempts; ++attempt)
	{
	  if (remaining_cmds.size() == 0)
	    break;
	  if (procs.size() == 0)
	  {
	    for (int i = 0; i < Layout::numNodes(); ++i)
	    {
	      procs.push_back(i);
	    }
	  }
	  std::vector<int> local_failed_cmds;
	  int proc = 0;
	  for (const auto& index : remaining_cmds)
	  {
	    if (procs.at(proc++ % procs.size()) == Layout::nodeNumber())
	    {
	      if (execute(params.param.cmd.at(index), params.param.constrain_to_gpu) != 0)
		local_failed_cmds.push_back(index);
	    }
	  }
#ifdef BUILD_REDSTAR_DATALIB
	  const auto& failed_cmds = SBN::gather(local_failed_cmds, SB::detail::getDefaultComm());
	  procs.resize(0);
	  for (int proc = 0; proc < failed_cmds.size(); ++proc)
	  {
	    if (failed_cmds.at(proc).size() == 0)
	    {
	      procs.push_back(proc);
	    }
	  }
	  remaining_cmds.resize(0);
	  for (const auto& indices : failed_cmds)
	  {
	    remaining_cmds.insert(remaining_cmds.end(), indices.begin(), indices.end());
	  }
#else
	  procs.resize(0);
	  procs.push_back(Layout::nodeNumber());
	  remaining_cmds = failed_cmds;
#endif
	}
	if (remaining_cmds.size() > 0)
	{
	  throw std::runtime_error("some tasks failed!");
	}
      }

      QMP_barrier();

      QDPIO::cout << InlineExecuteCmdEnv::name << ": ran successfully" << std::endl;

      END_CODE();
    }
  }
}
