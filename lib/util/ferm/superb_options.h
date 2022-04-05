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

      void get_prev_lines(const std::string& file, unsigned int char_num,
			  unsigned int num_prev_lines = 0, const std::string prefix = "") const
      {
	if (char_num > file.size())
	  throw std::runtime_error("Erroneous character number");
	std::size_t line_num = std::count(file.begin(), file.begin() + char_num, '\n') + 1;
	auto p = file.begin();
	for (unsigned int l = 1; l <= line_num && p != file.end(); ++l)
	{
	  auto p1 = std::find(p, file.end(), '\n');
	  if (l + num_prev_lines >= line_num)
	    QDPIO::cout << prefix << std::string(p, p1) << std::endl;
	  p = p1 + (p1 != file.end() ? 1 : 0);
	}
      }
    }

    /// Class for storing options
    struct Options {
      /// Get track of the path of this option on a set of options
      std::string prefix;
      /// Content of the file where the option comes from
      std::shared_ptr<std::string> file;
      /// First line number in `file` associated to the option
      unsigned int char_num;
      /// Whether this option has been checked
      mutable bool visited;

    protected:
      Options() : char_num(0), visited(false)
      {
      }

      /// Copy `prefix`, `file` and `char_num` from a given option
      /// \param op: option to copy the info

      void copyFileInfo(const Options& op)
      {
	  prefix = op.prefix;
	  file = op.file;
	  char_num = op.char_num;
      }

      /// Throw an error
      /// \param s: error message

      void throw_error(const std::string& err_msg) const
      {
	if (file)
	  throw std::runtime_error("Error at prefix `" + prefix + "'(l. " +
				   std::to_string(line_num) + "): " + err_msg + "\n" +
				   detail::get_prev_lines(*file, char_num, 5, "| "));
	else
	  throw std::runtime_error("Error at prefix `" + prefix + ": " + err_msg);
      }

    public:
      /// Type of the option
      enum Type { None, String, Double, Vector, Dictionary };

      /// Return the type of this option
      virtual Type getType() const = 0;

      /// Return if this options isn't None
      explicit operator bool() const noexcept
      {
	return getType() != None;
      }

      /// Return the string content of the option
      virtual std::string getString() const
      {
	throw_error("expected the value to be a string");
      }

      /// Return the double content of the option
      virtual double getDouble() const
      {
	throw_error("expected the value to be a double");
      }

      /// Return the integer content of the option
      virtual int getInt() const
      {
	throw_error("expected the value to be an integer");
      }

      /// Return the unsigned integer content of the option
      virtual unsigned int getUInt() const
      {
	throw_error("expected the value to be an unsigned integer");
      }

      /// Return the unsigned integer content of the option
      virtual bool getBool() const
      {
	throw_error("expected the value to be a boolean");
      }

      /// Return the vector content of the vector
      virtual const std::vector<Options>& getVector() const
      {
	throw_error("expected the value to be a vector");
      }

      /// Return the vector content of the vector
      virtual std::vector<Options>& getVector()
      {
	throw_error("expected the value to be a vector");
      }

      /// Return the map content of the vector
      virtual const std::map<std::string, Options>& getDictionary() const
      {
	throw_error("expected the value to be a dictionary");
      }

      /// Return the map content of the vector
      virtual std::map<std::string, Options>& getDictionary()
      {
	throw_error("expected the value to be a dictionary");
      }

      /// Return the option content on a path
      Option getValue(const std::string& path, Maybe<Option> defaultValue = none,
		      Maybe<Type> expectedType, Maybe<Options> fromOption = none,
		      std::string currentPrefix = "") const
      {
	const std::string errorHeader = "Error in searching for option `" + currentPrefix + path +
					"' from option at `" + fromOption.getSome(*this).prefix +
					"': ";

	if (path.size() == 0 || path == std::string("/")) {
	  if (expectedType && getType() != expectedType.getSome())
	    throw std::runtime_error(errorHeader + "Expected another type");
	  return *this;
	}
	if (getType() != Dictionary)
	{
	  if (defaultValue)
	    return defaultValue.getSome();
	  throw std::runtime_error(errorHeader + "the element `" + path "' is not a dictionary");
	}
	if (path[0] == '/')
	  return getValue(std::string(path.begin() + 1, path.end()), defaultValue, prefix + "/");
	auto p = std::find(path.begin(), path.end(), '/');
	std::string fieldName = std::string(path.begin(), p);
	auto m = getDictionary();
	if (m.count(fieldName) == 0) {
	  if (defaultValue.hasSome())
	    return defaultValue.getSome();
	  throw std::runtime_error(errorHeader + "the tag `" + fieldname + "' was not found");
	}
	return m[fieldName].getValue(std::string(p, path.end()), defaultValue,
				     prefix + std::string(path.begin(), p));
      }

      void setPrefix(const std::string& thisPrefix = "")
      {
	// Set the prefix of this option
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
	  usigned int i = 0;
	  for (auto& it : getVector())
	    it.setPrefix(thisPrefix + "[" + std::to_string(i++) + "]/");
	  break;
	}

	case Dictionary:
	{
	  for (auto& it : getDictionary())
	    it.second.setPrefix(thisPrefix + it.first + "/");
	  break;
	}
	}
      }
    };

    /// Storing a string as the value of an option
    struct NoneOption : public Options {
      NoneOption(const std::string& s, std::shared_ptr<std::string> file, unsigned int char_num)
	: file(file), char_num(char_num)
      {
      }
      Type getType() override const
      {
	return None;
      }
    };

    /// Storing a string as the value of an option
    struct StringOption : public Options {
      std::string value;
      StringOption(const std::string& s, std::shared_ptr<std::string> file, unsigned int char_num)
	: value(s), file(file), char_num(char_num)
      {
      }
      StringOption(const std::string& s, Maybe<Options> op=none) : value(s)
      {
	if (op)
	  copyFileInfo(op.getSome());
      }

      Type getType() override const
      {
	return String;
      }
      std::string getString() override const
      {
	visited = true;
	return value;
      }
      int getInt() override const
      {
	visited = true;
	try
	{
	  return std::stoi(value);
	} catch (...)
	{
	  throw_error("expected the value to be an integer");
	}
      }
      unsigned int getUInt() override const
      {
	visited = true;
	try
	{
	  return std::stoul(value);
	} catch (...)
	{
	  throw_error("expected the value to be an unsigned integer");
	}
      }
      bool getBool() override const
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
      }
      std::vector<Options> getVector() override const
      {
	visited = true;
	std::vector<Options> v;
	for (auto i = value.begin(), w = value.begin(); i != value.end(); ++i)
	{
	  if (std::isblank(*i) || i + 1 == value.end())
	  {
	    if (!std::isblank(*w))
	      v.push_back(StringOptions{std::string(w, i + 1 == value.end() ? i + 1 : i), *this});
	    w = i + 1;
	  }
	}
	return v;
      }
    };

    /// Storing a double as the value of an option
    struct DoubleOption : public Options {
      double value;
      DoubleOption(double d, std::shared_ptr<std::string> file, unsigned int char_num)
	: value(d), file(file), char_num(char_num)
      {
      }
      DoubleOption(double d, Maybe<Options> op=none) : value(d)
      {
	if (op)
	  copyFileInfo(op.getSome());
      }

      Type getType() override const
      {
	return Double;
      }
      double getDouble() override const
      {
	visited = true;
	return value;
      }
      int getInt() override const
      {
	visited = true;
	return (int)std::round(value);
      }
      unsigned int getUInt() override const
      {
	visited = true;
	int r = getInt();
	if (r < 0)
	  throw std::runtime_error("Error at prefix `" + prefix +
				   "': expected the value to be an unsigned integer");
	return (unsigned int)r;
      }
      bool getBool() override const
      {
	visited = true;
	return std::fabs(value) != 0;
      }
    };

    /// Storing a vector as the value of an option
    struct VectorOption : public Options {
      std::vector<Options> value;
      VectorOption(const std::vector<Options>& v, std::shared_ptr<std::string> file,
		   unsigned int char_num)
	: value(v), file(file), char_num(char_num)
      {
      }
      VectorOption(const std::vector<Options> &v, Maybe<Options> op=none) : value(v)
      {
	if (op)
	  copyFileInfo(op.getSome());
      }

      Type getType() override const
      {
	return Vector;
      }
      const std::vector<Options>& getVector() override const
      {
	visited = true;
	return value;
      }
      std::vector<Options>& getVector() override
      {
	visited = true;
	return value;
      }
    };

    /// Storing a vector as the value of an option
    struct DictionaryOption : public Options {
      std::map<std::string, Options> value;
      DictionaryOption(const std::map<std::string, Options>& m, std::shared_ptr<std::string> file,
		       unsigned int char_num)
	: value(m), file(file), char_num(char_num)
      {
      }
      DictionaryOption(const std::map<std::string, Options>& m, Maybe<Options> op=none) : value(m)
      {
	if (op)
	  copyFileInfo(op.getSome());
      }

      Type getType() override const
      {
	return Dictionary;
      }
      const std::map<std::string, Options>& getVector() override const
      {
	visited = true;
	return value;
      }
      std::map<std::string, Options>& getVector() override
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
    T getOption(const Options& ops, const std::string& path, Maybe<T> defaultValue = none);

    /// Return a string option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    std::string getOption<std::string>(const Options& ops, const std::string& path,
				       Maybe<std::string> defaultValue = none)
    {
      return ops.getValue(path, defaultValue ? StringOption{defaultValue.getSome()} : none)
	.getString();
    }

    /// Return a double option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    double getOption<double>(const Options& ops, const std::string& path,
			     Maybe<double> defaultValue = none)
    {
      return ops.getValue(path, defaultValue ? DoubleOption{defaultValue.getSome()} : none)
	.getDouble();
    }

    /// Return an integer option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    int getOption<int>(const Options& ops, const std::string& path, Maybe<int> defaultValue = none)
    {
      return ops.getValue(path, defaultValue ? DoubleOption{defaultValue.getSome()} : none)
	.getInt();
    }

    /// Return an unsigned integer option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    unsigned int getOption<unsigned int>(const Options& ops, const std::string& path,
			   Maybe<unsigned int> defaultValue = none)
    {
      return ops.getValue(path, defaultValue ? DoubleOption{defaultValue.getSome()} : none)
	.getUInt();
    }

    /// Return a boolean option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <>
    bool getOption<bool>(const Options& ops, const std::string& path,
			   Maybe<bool> defaultValue = none)
    {
      return ops
	.getValue(path,
		  defaultValue ? StringOption{defaultValue.getSome() ? "true" : "false"} : none)
	.getBool();
    }

    /// Return an enum option given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <typename Enum>
    Enum getOption(const Options& ops, const std::string& path,
		   const std::map<std::string, Enum>& m, Maybe<Enum> defaultValue = none)
    {
      const std::string defaultStr = "default";
      std::string value =
	ops.getValue(path, defaultValue ? StringOption{defaultStr} : none).getString();
      if (value == defaultStr)
	return defaultValue;
      if (m.count(value) == 0)
      {
	std::string availableOptions = "";
	for (const auto& it : m)
	  availableOptions += it.first + " ";
	ops.throw_error("unsupported value `" + value + "'; supported values: " + availableOptions);
      }
      return m[value];
    }

    /// Return a vector of options given a path
    /// \param ops: options into look for
    /// \param path: option path
    /// \param defaultValue: return value if the options isn't specified

    template <typename T>
    std::vector<T> getOption(const Options& ops, const std::string& path,
			     Maybe<std::vector<T>> defaultValue = none)
    {
      Option valueOp = ops.getValue(path, defaultValue ? NoneOption{} : none);
      if (!valueOp)
	return defaultValue;
      std::vector<T> r;
      for (const auto& op : valueOp.getVector())
	r.push_back(getOption<T>(op, ""));
      return r;
    }

    /// Returns options from XML
    /// \param s: text

    inline Options getOptionsFromXML(const std::string& s)
    {

      auto throw_error = [](std::string::iterator i, std::string err_msg) {
	std::size_t line_number = std::count(s.begin(), i, '\n') + 1;
	throw std::runtime_error("Error parsing XML at line " + std::to_string(line_number) + ": " +
				 err_msg);
      };

      std::vector<Options> ops;		// options of the open tags
      std::vector<std::string> tags;	// names of the open tags
      std::vector<Options::Type> types; // type of the open tags

      // Put the root element
      tags.push_back("/");
      types.push_back(Options::None);
      std::shared_ptr<std::string> file = std::make_shared(file);
      ops.push_back(NoneOptions{file, 0});

      // Parse the string
      for (auto i = s.begin(); i != s.end(); ++i)
      {
	// Ignore blanks
	if (std::isblank(*i))
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
	    auto i1 = std::find(i0, s.end(), '>');
	    if (i1 == s.end())
	      throw_error(i, "Broken tag: `<' without `>'");
	    std::string e("/>");
	    bool closing_tag = std::search(i0, i1, e.begin(), e.end()) != i1;

	    // Set the type of the current option
	    bool is_vector = tag == std::string("elem");
	    if (is_vector && types.last != None && types.last != Vector)
	      throw_error(i, "Mixing <elem> tags with other tags is not supported");
	    if (!is_vector && types.last != None && types.last != Dictionay)
	      throw_error(i, "Mixing <elem> tags with other tags is not supported");
	    if (types.last == Options::None)
	    {
	      if (is_vector)
		ops.last = VectorOption{{}, ops.last};
	      else
		ops.last = DictionaryOption{{}, ops.last};
	      types.last = (is_vector ? Options::Vector : Options::Dictionary);
	    }
	    else if (!is_vector && ops.last.getDictionary().count(tag) != 0)
	    {
	      throw_error(i,
			  "Unsupported more than one tag with the same name under the same node");
	    }

	    Options new_op = NoneOptions{file, i - s.begin()}; // new option to be filled
	    if (!closing_tag)
	    {
	      // Open tag
	      tags.push_back(tag);
	      types.push_back(NoneOption);
	      ops.push_back(new_op);
	    }
	    else
	    {
	      // Open and close tag
	      if (is_vector)
		ops.last.getVector().push_back(new_op);
	      else
		ops.last.getDictionary().insert({tag, new_op});
	    }
	  }
	  else if (*(i + 1) == '/')
	  {
	    // This is a close tag
	    auto i0 = i + 1;
	    for (; i0 != s.end() && std::isalpha(*i0); ++i0)
	      ;
	    std::string tag(i + 1, i0);
	    auto i1 = std::find(i0, s.end(), '>');
	    if (i1 == s.end())
	      throw_error(i, "Broken tag: `<' without `>'");
	    std::string e("/>");
	    bool closing_tag = std::search(i0, i1, e.begin(), e.end()) != i1;
	    if (closing_tag)
	      throw_error("Malformed tag `</' ... `/>'");
	    if (tag != tags.last)
	      throw_error(i, "Unmatched closing tag </" + tag + ">");

	    // Consume the option on top
	    types.pop_back();
	    Option op = ops.last;
	    ops.pop_back();
	    tags.pop_back();
	    if (type.last == Options::Vector)
	      ops.last.getVector().push_back(op);
	    else
	      ops.last.getDictionary().insert({tag, op});
	  }
	}

	// Detect content
	if (types.back() != Options::None)
	  throw_error(i, "Not supported text between tags");

	auto p = std::find(i, s.end(), '<'); // first non-content character
	auto i1 = p - 1;
	for (; std::isblank(*i1); --i1)
	  ;
	std::string content(i, i1 + 1);

	// Set content into the top option
	ops.last = StringOption{content, ops.last};
	i = p - 1;
      }

      if (tags.size() != 1)
	throw_error(s.end(), "there are unclosed tags");

      if (ops.size() == 0)
	return NoneOption{};

      return ops.first;
    }
  }
}

#endif // BUILD_SB

#endif // __INCLUDE_SUPERB_OPTIONS__
