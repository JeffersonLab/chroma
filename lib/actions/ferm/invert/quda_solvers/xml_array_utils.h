#ifndef XML_ARRAY_UTILS_H
#define XML_ARRAY_UTILS_H

namespace Chroma {
	// Read something optionally -- if not found leave it alone
	template<typename T>
	void readOptional(XMLReader& xml_in, const std::string& path, T& result)
	{
		if( xml_in.count(path) == 1 ) { 
			read(xml_in, path, result);	
		}
		else {
			QDPIO::cout << "Tag " << path << " not found\n";
		}
	}

	// Read something optionally -- if not found, set default valuej
	template<typename T>
	void readOptional(XMLReader& xml_in, const std::string& path, T& result, const T& defValue)
	{
		if ( xml_in.count(path) == 0 ) {
			QDPIO::cout << "Tag " << path << " not found using default: " << defValue << "\n";
			result = defValue;
		}
		else read(xml_in,path,result);
	}

	inline
	void readQudaBool(XMLReader& xml_in, const std::string&path, QudaBoolean& b)
	{
		bool tmp;
		read(xml_in, path, tmp);
		b = ( tmp == true ) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
	}

	inline
	void readQudaBool(XMLReader& xml_in, const std::string& path, QudaBoolean& b, bool defValue)
	{
		if( xml_in.count(path) == 0 ) b = ( defValue ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE );
		else readQudaBool(xml_in, path, b);
	}

	template<typename T>
	void readArray(XMLReader& paramtop, const std::string& path, multi1d<T>& array, const T& defValue)
	{

		multi1d<T> tmp;
		// If path is not found use default
		if ( paramtop.count(path) == 0 ) {

			QDPIO::cout << "Parameter with " << path << " not found. Setting default value "
					<< defValue <<  " for " << array.size() << " array members" << std::endl;

			for(int l=0; l < array.size() ; ++l) array[l] = defValue;
		}
		else {

			// If it is found read it to tmp
			read(paramtop, path, tmp);
			if ( tmp.size() == 1 ) {
                QDPIO::cout << "Broadcasting " << path << " = " << tmp[0] << "  to "  << array.size() <<  " array members" << std::endl;
				// if tmp is a single element array, broadcast it
				for(int l=0; l < array.size(); ++l) array[l] = tmp[0];
			}
			else {

				// If tmp is the same size as array copy it
				QDPIO::cout << "Copying " << path << " values to " << array.size() << " members " << std::endl;
				if ( tmp.size() == array.size() ) {
					for(int l=0; l < array.size(); ++l) array[l] = tmp[l];
				}
				else {
					QDPIO::cout << "Error: Array with path " << path << "has size "
							<< tmp.size() << " but " << array.size() << " are expected. " << std::endl;
					QDP_abort(1);
				}
			}
		}
	}

	inline
	void write(XMLWriter& xml, const std::string& path, QudaBoolean& b) {
		bool t = (b == QUDA_BOOLEAN_TRUE) ? true : false;
		write(xml, path, t);
	}
}

#endif