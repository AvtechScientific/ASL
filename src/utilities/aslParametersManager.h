/*
 * Advanced Simulation Library <http://asl.org.il>
 * 
 * Copyright 2015 Avtech Scientific <http://avtechscientific.com>
 *
 *
 * This file is part of Advanced Simulation Library (ASL).
 *
 * ASL is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, version 3 of the License.
 *
 * ASL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with ASL. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef ASLPARAMETERSMANAGER_H
#define ASLPARAMETERSMANAGER_H

#include "aslUValue.h"
#include <boost/program_options.hpp>
#include <map>
#include <typeinfo>

namespace asl
{

	class PrefixStore;

	/// This class stores parameter's value and the information
	/// needed to extract it from command line and/or parameters file
	/// \ingroup LDI
	template <typename T> class Parameter
	{
		public:
			/// \p key is the parameter's identification key in the
			/// parameters file. If no default value is specified, then the
			/// parameter is required to be specified in the parameters file
			/// or command line.
			Parameter(std::string key_,
			          std::string description_,
			          std::string units_ = "");
			/// \p key is the parameter's identification key in the parameters file
			/// If a default value is specified, then the parameter is not
			/// required to be specified in the parameters file or command line.	
			Parameter(T defaultValue,
			          std::string key_,
			          std::string description_,
			          std::string units_ = "");
			inline const T & v() const;
			inline T & v();
			inline std::shared_ptr<T> p();

		private:
			UValue<T> parameter;
			std::string key;
			std::string description;
			std::string units;
	};


	/** This class automatically accomodates newly created Parameters and then
	can load them from a parameters file.
	It has to be declared before declaring all the parameters it will manage!
	\ingroup Utilities */
	class ParametersManager
	{
		public:
			ParametersManager();
			~ParametersManager();
			/// Enables parameter loader
			void enable();
			/// Adds a Parameter to ParametersManager
			template <typename T> void add(UValue<T> parameter,
			                               std::string key,
			                               std::string description);
			/// Adds a group of parameters with common prefix to ParametersManager
			template <typename T> void add(UValue<std::map<std::string, T>> parameter,
			                               std::string key,
			                               std::string description);
			/// Adds a Parameter to ParametersManager
			template <typename T> void add(UValue<T> parameter,
			                               T defaultValue,
			                               std::string key,
			                               std::string description);
			/// Adds prefix and the pointer on the respective
			/// Parameter's destinationMap
			template <typename T>
			void addPrefix(const std::string prefix,
			               std::shared_ptr<std::map<std::string, T>> destinationMap);
			/// Loads all previously declared parameters
			/// from parameters file \p paramFile
			void load(std::string paramFile);
			std::string getFolder();
			std::string getFolderWithSlash();

			static ParametersManager * current;

		protected:
			boost::program_options::options_description parametersOptions;
			std::string folder;
			std::string folderWithSlash;
			/// Accomodates prefixes (defined by attached "*" wildcard)
			/// using PrefixStore class
			std::vector<std::shared_ptr<PrefixStore>> prefixes;

			void populateMaps(boost::program_options::variables_map & vm);
			/// Wrties all parameters and their
			/// default values (if available) to the file \p fileName
			void writeParametersFile(const std::string fileName);
			/// Content of the parameters file
			std::string parametersFileStr;
	};


	/** This class inherits ParametersManager class and thus also automatically
	accomodates newly created Parameters and then can load them from
	a parameters file and/or command line. It silently includes two parameters -
	`platform` and `device` that determine the hardware application will run on.
	It has to be declared before declaring all the parameters it will manage!
	\ingroup Utilities */
	class ApplicationParametersManager: public ParametersManager
	{
		public:
			ApplicationParametersManager(std::string applicationName_,
			                             std::string applicationVersion_,
			                             std::string paramFileName_ = "parameters.ini");
			
			/** Loads all previously declared parameters from command line
			and/or parameters file (provided through command line) */
			void load(int argc, char* argv[]);

		private:
			UValue<std::string> platform;
			UValue<std::string> device;
			std::string applicationName;
			std::string applicationVersion;
			std::string paramFileName;
			
	};


//-------------------------- Implementation --------------------------


	template <typename T> const T & Parameter<T>::v() const
	{
		return parameter.v();
	}


	template <typename T> T & Parameter<T>::v()
	{
		return parameter.v();
	}

	template <typename T> std::shared_ptr<T> Parameter<T>::p()
	{
		return parameter.p;
	}
	
} //namespace asl


#endif // ASLPARAMETERSMANAGER_H
