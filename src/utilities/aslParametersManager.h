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

	/** This class stores parameter's value and the information
	 needed to extract it from command line and/or parameters file.
	 Important: declare Parameters only after declaring ParametersManager or
	 ApplicationParametersManager instance and before calling 
	 ParametersManager::load() because each Parameter adds itself
	 to the instance automatically!
	\ingroup LDI */
	template <typename T> class Parameter
	{
		public:
			/**
			  \p key_ - option key; is used to specify this parameter through
			  command line and/or parameters file.
			  \p description_ is used in the help output and as comment on
			  parameters file generation.
			  \p units_ - parameter units; is used to complement the option
			  description mentioned above. Might be used for automatic unit
			  conversion in future (to this end it is recommended to use the
			  notation of the Boost::Units library).
			  Since no default value is specified, the parameter is required
			  to be specified in the parameters file or command line.
			*/
			Parameter(const char* key_,
			          const char* description_,
			          const char* units_ = "");
			/**
			  \p defaultValue - default parameter value. Since a default value is
			  specified, the parameter is not required to be specified in the
			  parameters file or command line.
			  \p key_ - option key; is used to specify this parameter through
			  command line and/or parameters file.
			  \p description_ is used in the help output and as comment on
			  parameters file generation.
			  \p units_ - parameter units; is used to complement the option
			  description mentioned above. Might be used for automatic unit
			  conversion in future (to this end it is recommended to use the
			  notation of the Boost::Units library).
			*/
			Parameter(T defaultValue,
			          const char* key_,
			          const char* description_,
			          const char* units_ = "");
			inline const T & v() const;
			inline T & v();
			inline std::shared_ptr<T> p();

		private:
			UValue<T> parameter;
			const std::string key;
			const std::string description;
			const std::string units;
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
			                               const char* key,
			                               const char* description,
			                               const char* units);
			/// Adds a group of parameters with common prefix to ParametersManager
			template <typename T> void add(UValue<std::map<std::string, T>> parameter,
			                               const char* key,
			                               const char* description,
			                               const char* units);
			/// Adds a Parameter with a default value to ParametersManager
			template <typename T> void add(UValue<T> parameter,
			                               T defaultValue,
			                               const char* key,
			                               const char* description,
			                               const char* units);
			/// Adds prefix and the pointer on the respective
			/// Parameter's destinationMap
			template <typename T>
			void addPrefix(const std::string prefix,
			               std::shared_ptr<std::map<std::string, T>> destinationMap);
			/// Loads all previously declared parameters
			/// from parameters file \p paramFile
			void load(std::string paramFile);
			/** Get the parameters file directory (absolute, with ending slash).
			It is convenient to write simulation results output together with
			its respective parameters file into the same directory (whose name
			might reflect the specifics of the parameters used) */
			std::string getDir();

			static ParametersManager * current;

		protected:
			boost::program_options::options_description parametersOptions;
			std::string parametersFileDirectory;
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
	`platform` and `device` that determine the hardware the application will run on.
	It has to be declared before declaring all the parameters it will manage!
	\ingroup Utilities */
	class ApplicationParametersManager: public ParametersManager
	{
		public:
			ApplicationParametersManager(const char* applicationName_,
			                             const char* applicationVersion_);
			
			/** Loads all previously declared parameters from command line
			and/or parameters file (provided through command line) */
			void load(int argc, char* argv[]);

		private:
			UValue<std::string> platform;
			UValue<std::string> device;
			std::string applicationName;
			std::string applicationVersion;
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
