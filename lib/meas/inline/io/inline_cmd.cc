/*! \file
 * \brief Inline task to execute commands
 */

#include "meas/inline/io/inline_cmd.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "qdp.h"
#include <cstdlib>

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
    }

    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "cmd", input.cmd);
      write(xml, "only_on_master", input.only_on_master);

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

    void execute(const std::string& cmd, int proc)
    {
      QDPIO::cout << "on proc " << proc << " executing " << cmd << std::endl;
      if (proc == Layout::nodeNumber())
      {
	int ret = std::system(cmd.c_str());
	if (ret == -1)
	  throw std::runtime_error("error executing a command line");
#if defined(_WIN32)
	const auto exit_code = ret;
#else
	const auto exit_code = WEXITSTATUS(ret);
#endif
	std::cout << "on proc " << proc << " result " << exit_code << std::endl;
      }
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
	for (const auto& cmd : params.param.cmd)
	  execute(cmd, 0);
      }
      else
      {
	if (params.param.cmd.size() > Layout::numNodes())
	  throw std::runtime_error("there are more commands that processes");
	int proc = 0;
	for (const auto& cmd : params.param.cmd)
	  execute(cmd, proc++);
      }

      QDPIO::cout << InlineExecuteCmdEnv::name << ": ran successfully" << std::endl;

      END_CODE();
    }
  }
}
