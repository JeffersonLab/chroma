// -*- C++ -*-
/*! \file                                                                    
 * \brief Multigrid prototype next
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#ifndef __INCLUDE_SUPERB_OPTIONS__
#define __INCLUDE_SUPERB_OPTIONS__

#include "chromabase.h"
#include "util/ferm/superb_contractions.h"

#ifdef BUILD_SB
namespace Chroma
{

  namespace SB
  {

    namespace detail {
      /// Return the previous lines up this one
      /// \param file: content
      /// \param char_num: character index of the last line to print
      /// \param num_prev_lines: maximum number of previous lines to print
      /// \param prefix: string to print previous to each line

      inline std::string get_prev_lines(const std::string& file, std::size_t char_num,
					unsigned int num_prev_lines = 0,
					const std::string prefix = "")
      {
	if (char_num > file.size())
	  throw std::runtime_error("Erroneous character number");
	std::size_t line_num = std::count(file.begin(), file.begin() + char_num, '\n') + 1;
	auto p = file.begin();
	std::string r;
	for (unsigned int l = 1; l <= line_num && p != file.end(); ++l)
	{
	  auto p1 = std::find(p, file.end(), '\n');
	  if (l + num_prev_lines >= line_num)
	    r += std::string(p, p1) + std::string("\n");
	  p = p1 + (p1 != file.end() ? 1 : 0);
	}
	return r;
      }
    }

    /// Class for storing options
    struct Option {
      /// Get track of the path of this option on a set of options
      std::string prefix;
      /// Content of the file where the option comes from
      std::shared_ptr<std::string> file;
      /// First line number in `file` associated to the option
      std::size_t char_num;
      /// Whether this option has been checked
      mutable bool visited;

      Option() : char_num(0), visited(false)
      {
      }

    protected:
      Option(std::shared_ptr<std::string> file, std::size_t char_num)
	: file(file), char_num(char_num), visited(false)
      {
      }

      /// Copy `prefix`, `file` and `char_num` from a given option
      /// \param op: option to copy the info

      void copyFileInfo(const Option& op)
      {
	  prefix = op.prefix;
	  file = op.file;
	  char_num = op.char_num;
      }

    public:

      /// Throw an error
      /// \param s: error message

      void throw_error(const std::string& err_msg) const
      {
	if (file)
	{
	  std::size_t line_num = std::count(file->begin(), file->begin() + char_num, '\n') + 1;
	  throw std::runtime_error(std::string("Error at prefix `") + prefix +	      //
				   "'(l. " + std::to_string(line_num) + "): " +	      //
				   err_msg + "\n" +				      //
				   detail::get_prev_lines(*file, char_num, 5, "| ")); //
	}
	else
	{
	  throw std::runtime_error("Error at prefix `" + prefix + ": " + err_msg);
	}
      }

      /// Type of the option
      enum Type { None, String, Double, Vector, Dictionary };

      /// Return the type of this option
      virtual Type getType() const
      {
        throw_error("getType: invalid object, it's abstract");
        throw std::exception{}; // silent no return warning
      }

      /// Return if this options isn't None
      explicit operator bool() const noexcept
      {
	return getType() != None;
      }

      /// Return the string content of the option
      virtual std::string getString() const
      {
	throw_error("expected the value to be a string");
	throw std::exception{}; // silent no return warning
      }

      /// Return the double content of the option
      virtual double getDouble() const
      {
	throw_error("expected the value to be a double");
	throw std::exception{}; // silent no return warning
      }

      /// Return the integer content of the option
      virtual int getInt() const
      {
	throw_error("expected the value to be an integer");
	throw std::exception{}; // silent no return warning
      }

      /// Return the unsigned integer content of the option
      virtual unsigned int getUInt() const
      {
	throw_error("expected the value to be an unsigned integer");
	throw std::exception{}; // silent no return warning
      }

      /// Return the unsigned integer content of the option
      virtual bool getBool() const
      {
	throw_error("expected the value to be a boolean");
	throw std::exception{}; // silent no return warning
      }

      /// Return the vector content of the vector
      virtual std::vector<std::shared_ptr<Option>> getVector() const
      {
	throw_error("expected the value to be a vector");
	throw std::exception{}; // silent no return warning
      }

      /// Return the vector content of the vector
      virtual std::vector<std::shared_ptr<Option>>& getVector()
      {
	throw_error("expected the value to be a vector");
	throw std::exception{}; // silent no return warning
      }

      /// Return the map content of the vector
      virtual const std::map<std::string, std::shared_ptr<Option>>& getDictionary() const
      {
	throw_error("expected the value to be a dictionary");
	throw std::exception{}; // silent no return warning
      }

      /// Return the map content of the vector
      virtual std::map<std::string, std::shared_ptr<Option>>& getDictionary()
      {
	throw_error("expected the value to be a dictionary");
	throw std::exception{}; // silent no return warning
      }

      /// Return the option content on a path
      const Option& getValue(const std::string& path, Maybe<const Option&> defaultValue = none,
			     Maybe<Type> expectedType = none,
			     Maybe<const Option&> fromOption = none,
			     Maybe<std::string> originalPath = none) const
      {
	// If fromOption is none, set this node
	Maybe<const Option&> fromOption_{fromOption.getSome(*this)};

	// If originalPath is none, set the given path
	Maybe<std::string> originalPath_{originalPath.getSome(path)};

	// Construct a nice error message
	const std::string errorHeader = "Error in searching for option `" +
					originalPath_.getSome() + "' from option at `" +
					fromOption_.getSome().prefix + "': ";

	// If the path is empty or ask for the root node, just return this node
	if (path.size() == 0 || path == std::string("/")) {
	  if (expectedType && getType() != expectedType.getSome())
	    throw std::runtime_error(errorHeader + "Expected another type");
	  return *this;
	}

	// If that path has a tag but the current option isn't a dictionary, either return
	// the default value or throw an error
	if (getType() != Dictionary)
	{
	  if (defaultValue)
	    return defaultValue.getSome();
	  throw std::runtime_error(errorHeader + "the element `" + path + "' is not a dictionary");
	}

	// If path starts with `/`, consume it and continue
	if (path[0] == '/')
	  return getValue(std::string(path.begin() + 1, path.end()), defaultValue, expectedType,
			  fromOption_, originalPath_);

	// Find the name of the tag
	auto p = std::find(path.begin(), path.end(), '/');
	std::string fieldName = std::string(path.begin(), p);

	// If the tag isn't under the current node, either return the default value or throw an error
	auto m = getDictionary();
	if (m.count(fieldName) == 0) {
	  if (defaultValue.hasSome())
	    return defaultValue.getSome();
	  throw std::runtime_error(errorHeader + "the tag `" + fieldName + "' was not found");
	}

	// Otherwise, consume the tag name, and continue
	return m[fieldName]->getValue(std::string(p, path.end()), defaultValue, expectedType,
				      fromOption_, originalPath_);
      }

      /// Return the option content on a path
      const Option& getValue(const std::string& path, Type expectedType) const
      {
	return getValue(path, none, expectedType);
      }

      /// Return the option content on a path
      Maybe<const Option&> getValueMaybe(const std::string& path, Maybe<Type> expectedType = none) const
      {
	struct AuxNone : Option {
	  Type getType() const override
	  {
	    return None;
	  }
	} defaultOp;
	const Option& op = getValue(path, Maybe<const Option&>{defaultOp}, expectedType);
	if (&op == &defaultOp)
	  return none;
	return op;
      }

      void setPrefix(const std::string& thisPrefix = "")
      { // Set the prefix of this option
	prefix = thisPrefix;

	switch (getType())
	{
	case None:
	case String:
	case Double:
	{
	  // Do nothing
	  break;
	}

	case Vector:
	{
	  unsigned int i = 0;
	  for (auto& it : getVector())
	    it->setPrefix(thisPrefix + "[" + std::to_string(i++) + "]/");
	  break;
	}

	case Dictionary:
	{
	  for (auto& it : getDictionary())
	    it.second->setPrefix(thisPrefix + it.first + "/");
	  break;
	}
	}
      }
    };

    /// Storing a string as the value of an option
    struct NoneOption : public Option {
      NoneOption()
      {
      }
      NoneOption(std::shared_ptr<std::string> file, std::size_t char_num) : Option{file, char_num}
      {
      }
      Type getType() const override
      {
	return None;
      }
    };

    /// Storing a string as the value of an option
    struct StringOption : public Option {
      std::string value;
      StringOption(const std::string& s, std::shared_ptr<std::string> file, std::size_t char_num)
	: Option{file, char_num}, value(s)
      {
      }
      StringOption(const std::string& s, Maybe<const Option&> op=none) : value(s)
      {
	if (op)
	  copyFileInfo(op.getSome());
      }

      Type getType() const override
      {
	return String;
      }
      std::string getString() const override
      {
	visited = true;
	return value;
      }
      int getInt() const override
      {
	visited = true;
	try
	{
	  return std::stoi(value);
	} catch (...)
	{
	  throw_error("expected the value to be an integer");
	}
	throw std::exception{}; // silent no return warning
      }
      unsigned int getUInt() const override
      {
	visited = true;
	try
	{
	  return std::stoul(value);
	} catch (...)
	{
	  throw_error("expected the value to be an unsigned integer");
	}
	throw std::exception{}; // silent no return warning
      }
      bool getBool() const override
      {
	visited = true;
	std::string valueLower = value;
	std::transform(valueLower.begin(), valueLower.end(), valueLower.begin(),
		       [](unsigned char c) { return std::tolower(c); });
	if (valueLower == std::string("true"))
	  return true;
	if (valueLower == std::string("false"))
	  return false;
	throw_error("expected the value to be boolean, either `true' or `false'");
	throw std::exception{}; // silent no return warning
      }
      std::vector<std::shared_ptr<Option>> getVector() const override
      {
	visited = true;
	std::vector<std::shared_ptr<Option>> v;
	for (auto i = value.begin(), w = value.begin(); i != value.end(); ++i)
	{
	  if (std::isspace(*i) || i + 1 == value.end())
	  {
	    if (!std::isspace(*w))
	      v.push_back(std::make_shared<StringOption>(
		std::string(w, i + 1 == value.end() ? i + 1 : i), *this));
	    w = i + 1;
	  }
	}
	return v;
      }
    };

    /// Storing a double as the value of an option
    struct DoubleOption : public Option {
      double value;
      DoubleOption(double d, std::shared_ptr<std::string> file, std::size_t char_num)
	: Option{file, char_num}, value(d)
      {
      }
      DoubleOption(double d, Maybe<const Option&> op=none) : value(d)
      {
	if (op)
	  copyFileInfo(op.getSome());
      }

      Type getType() const override
      {
	return Double;
      }
      double getDouble() const override
      {
	visited = true;
	return value;
      }
      int getInt() const override
      {
	visited = true;
	return (int)std::round(value);
      }
      unsigned int getUInt() const override
      {
	visited = true;
	int r = getInt();
	if (r < 0)
	  throw std::runtime_error("Error at prefix `" + prefix +
				   "': expected the value to be an unsigned integer");
	return (unsigned int)r;
      }
      bool getBool() const override
      {
	visited = true;
	return std::fabs(value) != 0;
      }
    };

    /// Storing a vector as the value of an option
    struct VectorOption : public Option {
      std::vector<std::shared_ptr<Option>> value;
      VectorOption(const Option& op)
      {
	copyFileInfo(op);
      }

      Type getType() const override
      {
	return Vector;
      }
      std::vector<std::shared_ptr<Option>> getVector() const override
      {
	visited = true;
	return value;
      }
      std::vector<std::shared_ptr<Option>>& getVector() override
      {
	visited = true;
	return value;
      }
    };

    /// Storing a vector as the value of an option
    struct DictionaryOption : public Option {
      std::map<std::string, std::shared_ptr<Option>> value;
      DictionaryOption(const Option& op) 
      {
	copyFileInfo(op);
      }

      Type getType() const override
      {
	return Dictionary;
      }
      const std::map<std::string, std::shared_ptr<Option>>& getDictionary() const override
      {
	visited = true;
	return value;
      }
      std::map<std::string, std::shared_ptr<Option>>& getDictionary() override
      {
	visited = true;
	return value;
      }
    };

    /// Return an option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <typename T>
    T getOption(const Option& ops, const std::string& path, Maybe<T> defaultValue = none);

    /// Return a string option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    inline std::string getOption<std::string>(const Option& ops, const std::string& path,
					      Maybe<std::string> defaultValue)
    {
      StringOption defaultOp{defaultValue.getSome("")};
      return ops.getValue(path, defaultValue ? Maybe<const Option&>{defaultOp} : none).getString();
    }

    /// Return a double option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    inline double getOption<double>(const Option& ops, const std::string& path,
				    Maybe<double> defaultValue)
    {
      DoubleOption defaultOp{defaultValue.getSome(0.0)};
      return ops.getValue(path, defaultValue ? Maybe<const Option&>{defaultOp} : none).getDouble();
    }

    /// Return an integer option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    inline int getOption<int>(const Option& ops, const std::string& path, Maybe<int> defaultValue)
    {
      DoubleOption defaultOp{(double)defaultValue.getSome(0)};
      return ops.getValue(path, defaultValue ? Maybe<const Option&>{defaultOp} : none).getInt();
    }

    /// Return an unsigned integer option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    inline unsigned int getOption<unsigned int>(const Option& ops, const std::string& path,
						Maybe<unsigned int> defaultValue)
    {
      DoubleOption defaultOp{(double)defaultValue.getSome(0)};
      return ops.getValue(path, defaultValue ? Maybe<const Option&>{defaultOp} : none).getUInt();
    }

    /// Return a boolean option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    inline bool getOption<bool>(const Option& ops, const std::string& path,
				Maybe<bool> defaultValue)
    {
      StringOption defaultOp{defaultValue.getSome(false) ? "true" : "false"};
      return ops.getValue(path, defaultValue ? Maybe<const Option&>{defaultOp} : none).getBool();
    }

    /// Return an enum option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <typename Enum>
    Enum getOption(const Option& ops, const std::string& path,
		   const std::map<std::string, Enum>& m, Maybe<Enum> defaultValue)
    {
      const std::string defaultStr = "default";
      StringOption defaultOp{defaultStr};
      std::string value =
	ops.getValue(path, defaultValue ? Maybe<const Option&>{defaultOp} : none).getString();
      if (value == defaultStr)
	return defaultValue.getSome();
      if (m.count(value) == 0)
      {
	std::string availableOption = "";
	for (const auto& it : m)
	  availableOption += it.first + " ";
	ops.throw_error("unsupported value `" + value + "'; supported values: " + availableOption);
      }
      return m.at(value);
    }

    /// Return a vector of options given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <typename T>
    std::vector<T> getVectorOption(const Option& ops, const std::string& path,
				   Maybe<std::vector<T>> defaultValue = none)
    {
      NoneOption defaultOp{};
      const Option& valueOp = ops.getValue(path, defaultValue ? Maybe<const Option&>{defaultOp} : none);
      if (!valueOp)
	return defaultValue.getSome();
      std::vector<T> r;
      for (const auto& op : valueOp.getVector())
	r.push_back(getOption<T>(*op, ""));
      return r;
    }

    /// Return a dictionary option
    /// \param ops: options into look for
    /// \param path: option path

    inline const Option& getOptions(const Option& ops, const std::string& path)
    {
      return ops.getValue(path, Option::Dictionary);
    }

    /// Return a dictionary option
    /// \param ops: options into look for
    /// \param path: option path

    inline Maybe<const Option&> getOptionsMaybe(const Option& ops, const std::string& path)
    {
      return ops.getValueMaybe(path, Option::Dictionary);
    }

    using Options = Option;

    /// Returns options from XML
    /// \param s: text

    inline std::shared_ptr<Option> getOptionsFromXML(const std::string& s)
    {

      auto throw_error = [&](std::string::const_iterator i, std::string err_msg) {
	std::size_t line_number = std::count(s.begin(), i, '\n') + 1;
	throw std::runtime_error("Error parsing XML at line " + std::to_string(line_number) + ": " +
				 err_msg);
      };

      std::vector<std::shared_ptr<Option>> ops;		// options of the open tags
      std::vector<std::string> tags;	// names of the open tags
      std::vector<Option::Type> types; // type of the open tags

      // Put the root element
      tags.push_back("/");
      types.push_back(Option::None);
      std::shared_ptr<std::string> file = std::make_shared<std::string>(s);
      ops.push_back(std::make_shared<NoneOption>(file, 0));

      // Parse the string
      for (auto i = s.begin(); i != s.end(); ++i)
      {
	// Ignore blanks
	if (std::isspace(*i))
	  continue;

	// Detect a token
	if (*i == '<')
	{
	  if (i + 1 == s.end())
	  {
	    throw_error(i, "Broken tag: `<' without closing `>'");
	  }
	  else if (*(i + 1) == '?')
	  {
	    // Ignore <? ... ?>
	    const std::string e("?>");
	    auto p = std::search(i + 2, s.end(), e.begin(), e.end());
	    if (p == s.end())
	      throw_error(i, "Broken tag: `<?' without closing `?>'");
	    i = p + (e.size() - 1);
	  }
	  else if (*(i + 1) == '!' && i + 2 != s.end() && i + 3 != s.end() && *(i + 2) == '-' &&
		   *(i + 3) == '-')
	  {
	    // Ignore <!-- ... -->
	    const std::string e("-->");
	    auto p = std::search(i, s.end(), e.begin(), e.end());
	    if (p == s.end())
	      throw_error(i, "Broken tag: `<!--' without closing `-->'");
	    i = p + (e.size() - 1);
	  }
	  else if (*(i + 1) != '/')
	  {
	    // This is an open tag
	    auto i0 = i + 1;
	    for (; i0 != s.end() && std::isalpha(*i0); ++i0)
	      ;
	    std::string tag(i + 1, i0);
	    if (tag.size() == 0)
	      throw_error(i, "Error: empty tag");
	    auto i1 = std::find(i0, s.end(), '>');
	    if (i1 == s.end())
	      throw_error(i, "Broken tag: `<' without `>'");
	    std::string e("/>");
	    bool closing_tag = std::search(i0, i1, e.begin(), e.end()) != i1;

	    // Set the type of the current option
	    bool is_vector = tag == std::string("elem");
	    if (is_vector && types.back() != Option::None && types.back() != Option::Vector)
	      throw_error(i, "Mixing <elem> tags with other tags is not supported");
	    if (!is_vector && types.back() != Option::None && types.back() != Option::Dictionary)
	      throw_error(i, "Mixing <elem> tags with other tags is not supported");
	    if (types.back() == Option::None)
	    {
	      if (is_vector)
		ops.back() = std::make_shared<VectorOption>(*ops.back());
	      else
		ops.back() = std::make_shared<DictionaryOption>(*ops.back());
	      types.back() = (is_vector ? Option::Vector : Option::Dictionary);
	    }
	    else if (!is_vector && ops.back()->getDictionary().count(tag) != 0)
	    {
	      throw_error(i,
			  "Unsupported more than one tag with the same name under the same node");
	    }

	    std::shared_ptr<Option> new_op = std::make_shared<NoneOption>(
	      file, (std::size_t)(i - s.begin())); // new option to be filled
	    if (!closing_tag)
	    {
	      // Open tag
	      tags.push_back(tag);
	      types.push_back(Option::None);
	      ops.push_back(new_op);
	    }
	    else
	    {
	      // Open and close tag
	      if (is_vector)
		ops.back()->getVector().push_back(new_op);
	      else
		ops.back()->getDictionary().insert({tag, new_op});
	    }
	    i = i1;
	  }
	  else if (*(i + 1) == '/')
	  {
	    // This is a close tag
	    auto i0 = i + 2;
	    for (; i0 != s.end() && std::isalpha(*i0); ++i0)
	      ;
	    std::string tag(i + 2, i0);
	    if (tag.size() == 0)
	      throw_error(i, "Error: empty tag");
	    auto i1 = std::find(i0, s.end(), '>');
	    if (i1 == s.end())
	      throw_error(i, "Broken tag: `<' without `>'");
	    std::string e("/>");
	    bool closing_tag = std::search(i0, i1, e.begin(), e.end()) != i1;
	    if (closing_tag)
	      throw_error(i, "Malformed tag `</' ... `/>'");
	    if (tag != tags.back())
	      throw_error(i, "Unmatched closing tag </" + tag + ">");

	    // Consume the option on top
	    types.pop_back();
	    std::shared_ptr<Option> op = ops.back();
	    ops.pop_back();
	    tags.pop_back();
	    if (types.back() == Option::Vector)
	      ops.back()->getVector().push_back(op);
	    else
	      ops.back()->getDictionary().insert({tag, op});
	    i = i1;
	  }
	}
	else
	{
	  // Detect content
	  if (types.back() != Option::None)
	    throw_error(i, "Not supported text between tags");

	  auto p = std::find(i, s.end(), '<'); // first non-content character
	  auto i1 = p - 1;
	  for (; std::isspace(*i1); --i1)
	    ;
	  std::string content(i, i1 + 1);

	  // Set content into the top option
	  ops.back() = std::make_shared<StringOption>(content, *ops.back());
	  i = p - 1;
	}
      }

      if (tags.size() != 1)
	throw_error(s.end(), "there are unclosed tags");

      if (ops.size() == 0)
	return std::make_shared<NoneOption>();

      return ops.front();
    }
  }
}

#endif // BUILD_SB

#endif // __INCLUDE_SUPERB_OPTIONS__
