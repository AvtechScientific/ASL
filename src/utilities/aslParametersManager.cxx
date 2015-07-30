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

#include "aslParametersManager.h"
#include "../aslUtilities.h"
#include "../acl/aclHardware.h"
#include "math/aslVectorsDynamicLength.h"
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include "acl/aclTypesList.h"
#include <map>
#include <memory>

using namespace std;
using namespace boost::filesystem;
using namespace boost::program_options;

namespace asl
{

	// Auxiliary infrastructure in order to enable
	// usage of "*" wildcard in options key and
	// map population.
	// This should be removed once program_options
	// enables population of a map<string, T>.
	//
	//
	// see also (add after the implementation bellow):
	// http://stackoverflow.com/questions/15554842/boostprogram-options-parameters-with-a-fixed-and-a-variable-token
	//
	class PrefixStore
	{
		public:
			PrefixStore(const string prefix_);
			virtual void store(variables_map & vm) = 0;

		protected:
			string prefix;
	};


	template <typename T> class PrefixMapStore : public PrefixStore
	{
		public:
			PrefixMapStore(const string prefix_,
			               shared_ptr<map<string, T>> destinationMap_);
			virtual void store(variables_map & vm);

		private:
			shared_ptr<map<string, T>> destinationMap;
	};


	PrefixStore::PrefixStore(const string prefix_):
		prefix(prefix_)
	{
	}


	template <typename T>
	PrefixMapStore<T>::PrefixMapStore(const string prefix_,
	                                  shared_ptr<map<string, T>> destinationMap_ ):
		PrefixStore(prefix_),
		destinationMap(destinationMap_)
	{
	}


	template <typename T>
	void PrefixMapStore<T>::store(variables_map & vm)
	{
		variables_map::iterator it;
		for (it = vm.begin(); it != vm.end(); ++it)
		{
			if ((*it).first.find(prefix) != string::npos)
			{
				pair<string, T> p(it->first, it->second.as<T>());
				destinationMap->insert(p);
			}
		}
	}
	// End of auxiliary infrastructure.


	template <typename T> class WildcardCheck
	{
		public:
			static void check(string key, shared_ptr<T> dummyPointer)
			{
				if (*key.rbegin() == '*')
					errorMessage("Parameter<T>::Parameter() - attempt to use \"*\" wildcard in the option key without providing corresponding map");
			}
	};


	template <typename T> class WildcardCheck<map<string, T>>
	{
		public:
			static void check(string key, shared_ptr<map<string, T>> destinationMap)
			{
				if (*key.rbegin() != '*')
					errorMessage("Parameter<map<string, T>>::Parameter() - no \"*\" wildcard in the option key");		

				// generate prefix by cutting option key's last char - "*"
				ParametersManager::current->addPrefix(key.substr(0, key.size() - 1),
				                                      destinationMap);
			}
	};
		

	template <typename T>
	Parameter<T>::Parameter(const char* key_,
	                        const char* description_,
	                        const char* units_):
		key(key_),
		description(description_),
		units(units_)
	{
		// Adds itself to current ParametersManager
		if (ParametersManager::current == NULL)
			errorMessage("ParametersManager was not instantiated and is not available");

		WildcardCheck<T>::check(key, parameter.p);

		ParametersManager::current->add(parameter,
		                                key.c_str(),
		                                description.c_str(),
		                                units.c_str());
	}

	template Parameter<string>::Parameter(const char* key,
	                                      const char* description,
	                                      const char* units);
		
	#define BOOST_TT_rep_expression(r, data, t) \
		template Parameter<t>::Parameter(const char* key, \
										 const char* description, \
										 const char* units); \
		template Parameter<AVec<t>>::Parameter(const char* key, \
											   const char* description, \
											   const char* units); \
		template Parameter<map<string, t>>::Parameter(const char* key, \
													  const char* description, \
													  const char* units); \
		template Parameter<map<string, AVec<t>>>::Parameter(const char* key, \
															const char* description, \
															const char* units);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression


	template <typename T> Parameter<T>::Parameter(T defaultValue,
	                                              const char* key_,
	                                              const char* description_,
	                                              const char* units_):
		key(key_),
		description(description_),
		units(units_)
	{
		// Adds itself to current ParametersManager
		if (ParametersManager::current == NULL)
			errorMessage("ParametersManager was not instantiated and is not available");
		
		ParametersManager::current->add(parameter,
		                                defaultValue,
		                                key.c_str(),
		                                description.c_str(),
		                                units.c_str());
	}

	template Parameter<string>::Parameter(string defaultValue,
	                                      const char* key,
	                                      const char* description,
	                                      const char* units);
		
	#define BOOST_TT_rep_expression(r, data, t) \
		template Parameter<t>::Parameter(t defaultValue, \
										 const char* key, \
										 const char* description, \
										 const char* units); \
		template Parameter<AVec<t>>::Parameter(AVec<t> defaultValue, \
											   const char* key, \
											   const char* description, \
											   const char* units);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression


	ParametersManager * ParametersManager::current(NULL);

	ParametersManager::ParametersManager():
		parametersOptions("Parameters options"),
		parametersFileStr("# Parameters file\n"),
		parametersFileDirectory("")
	{
		enable();
	}


	ParametersManager::~ParametersManager()
	{
		// Deactivates this instance of ParametersManager
		if (current == this)
			current = NULL;
	}


	void ParametersManager::enable()
	{
		// Activates this instance of ParametersManager
		current = this;
	}


	// Function that teaches boost::program_options how to deal with AVec<>
	template <typename T> void validate(boost::any& v,
	                                    const std::vector<std::string>& values,
	                                    AVec<T> *, int)
	{
		// Make sure no previous assignment to 'v' was made.
		validators::check_first_occurrence(v);		

		AVec<T> vec;
		// Extract tokens from values string vector and populate AVec
		if (values[0] == "")
		{
			errorMessage("ParametersManager - no value provided for a variable of type AVec<...>");
		}

		vector<T> converted;
		stringstream strStream(values[0]);
		string token;
		while (!strStream.eof())
		{
			strStream >> token;
			converted.push_back(strToNum<T>(token));
		}

		vec = converted;
		v = vec;
	}


	template <typename T> void ParametersManager::add(UValue<T> parameter,
	                                                  const char* key,
	                                                  const char* description,
	                                                  const char* units)
	{
		string fullDescription = units[0] == '\0' ? description
								: string(description) + ", [" + units + "]";
		parametersOptions.add_options()
			(key, value<T>(&parameter.v())->required(), fullDescription.c_str());

		// Add the parameter to the default parameters file			
		parametersFileStr += "\n# " + fullDescription + "\n" + key + " = \n";
	}


	template<typename T> void ParametersManager::add(UValue<map<string, T>> parameter,
	                                                 const char* key,
	                                                 const char* description,
	                                                 const char* units)
	{
		string fullDescription = units[0] == '\0' ? description
								: string(description) + ", [" + units + "]";
		parametersOptions.add_options()
			(key, value<T>()->required(), fullDescription.c_str());
		// Add the parameter to the default parameters file
		parametersFileStr += "\n# " + fullDescription + "; replace '*' by any suffix to provide a set of related parameters\n"
							+ key + " = \n";
	}


	template<typename T> void ParametersManager::add(UValue<T> parameter,
	                                                 T defaultValue,
	                                                 const char* key,
	                                                 const char* description,
	                                                 const char* units)
	{
		string fullDescription = units[0] == '\0' ? description
								: string(description) + ", [" + units + "]";
		parametersOptions.add_options()
			(key, value<T>(&parameter.v())->default_value(defaultValue),
			 fullDescription.c_str());

		// Add the parameter to the default parameters file
		parametersFileStr += "\n# " + fullDescription + "\n"
							+ key + " = " + numToStr(defaultValue) + "\n";
	}


	template<typename T>
	void ParametersManager::addPrefix(string prefix,
	                                  shared_ptr<map<string, T>> destinationMap)
	{
		prefixes.push_back(std::make_shared<PrefixMapStore<T>>(prefix, destinationMap));
	}


	void ParametersManager::populateMaps(variables_map & vm)
	{
		for (unsigned int i = 0; i < prefixes.size(); ++i)
		{
			prefixes[i]->store(vm);
		}
	}


	void ParametersManager::load(string paramFile)
	{
		variables_map vm;

		try
		{
			ifstream ifs(paramFile);
			if (!ifs.good())
				errorMessage("Can not open parameters file: " + paramFile);

			parsed_options parsed = parse_config_file(ifs, parametersOptions,
			                                          true);
			store(parsed, vm);
			notify(vm);
			populateMaps(vm);
	}
		catch(exception& e)
		{
			errorMessage(string("ParametersManager::load() - ") + e.what());
		}
	}


	string ParametersManager::getDir()
	{
		return parametersFileDirectory;
	}
	

	void ParametersManager::writeParametersFile(const std::string fileName)
	{
		ofstream fo(fileName);
		if (!fo.good())
			errorMessage("ParametersManager::writeParametersFile() - can not open file: " + fileName);

		fo << parametersFileStr;
		fo.close();
	}


	ApplicationParametersManager::ApplicationParametersManager(const char* applicationName_,
	                                                           const char* applicationVersion_):
		platform(acl::getPlatformVendor(acl::hardware.defaultQueue)),
		device(acl::getDeviceName(acl::hardware.defaultQueue)),
		applicationName(applicationName_),
		applicationVersion(applicationVersion_)
	{
		enable();
		// Prepend informative header
		parametersFileStr += "# Generated by " + applicationName + " version " + applicationVersion + "\n\n"
							+ "# Get the list of all available computation devices by running:\n"
							+ "# `" + applicationName + " -d`\n";

		/* Add platform and device parameters already initialized
		with their default values */
		add(platform, platform.v(),
		    "platform", "Default computation platform", "");
		add(device, device.v(),
		    "device", "Default computation device", "");
	}


	void ApplicationParametersManager::load(int argc, char * argv[])
	{
		variables_map vm;
		// Set default parameters file path
		path p("./");

		options_description genericOptions("Generic options");

		genericOptions.add_options()
			("help,h", "Display this help and exit")
			("version,v", "Display version and exit")
			("devices,d", "Display all available devices and exit")
			("parameters,p", value<string>(), "Path to the parameters file")
			("generate,g", value<string>(),
			 "Generate default parameters file, write it and exit")
			("check,c", "Check parameters for consistency and exit");

		positional_options_description positional;

		options_description allOptions;
		positional.add("parameters", 1);

		allOptions.add(genericOptions).add(parametersOptions);

		try
		{
			store(command_line_parser(argc, argv).options(allOptions).positional(positional).run(), vm);

			if (vm.count("help"))
			{
				cout << "Usage: " + applicationName + " [PARAMETERS_FILE] [OPTION]...\n"
					 << allOptions
					 << endl;
				exit(EXIT_SUCCESS);
			}

			if (vm.count("version"))
			{
				cout << applicationName + " " + applicationVersion
					 << endl;
				exit(EXIT_SUCCESS);
			}

			if (vm.count("devices"))
			{
				cout << applicationName + " " + applicationVersion + "\n\n"
					<< "Default computation device:\n"
					<< acl::hardware.getDefaultDeviceInfo() << "\n\n"
					<< "List of all available platforms and their devices:\n"
					<< acl::hardware.getDevicesInfo()
					<< endl;
				exit(EXIT_SUCCESS);
			}

			if (vm.count("generate"))
			{
				path gp(vm["generate"].as<string>());
				cout << "Writing default parameters file to: "
					 << gp.string() << endl;

				writeParametersFile(gp.string());
				exit(EXIT_SUCCESS);
			}

			if (vm.count("parameters"))
			{
				p = vm["parameters"].as<string>();
				if (is_directory(p))
				{
					// Only warn, since all options might have default values, or required values
					// provided through the command line - so no parameters file is required
					warningMessage("ApplicationParametersManager::load() - no parameters file provided, " + p.string() + " is a directory. Using default and/or command line values");
				}
				else
				{
					ifstream ifs(p.string());
					if (ifs.good())
					{
						parsed_options parsed = parse_config_file(ifs, allOptions, true);
						store(parsed, vm);
						// Get directory, cutting file name
						p = p.parent_path();
					}
					else
					{
						warningMessage("ApplicationParametersManager::load() - can not open parameters file: " + p.string() + " . Using default and/or command line values");
						// Get directory, cutting file name, after using p for warning
						p = p.parent_path();
					}
				}
			}
			else
			{
				warningMessage("ApplicationParametersManager::load() - no parameters file provided. Using default and/or command line values");
			}

			// Generate `parametersFileDirectory`
			// Get absolute, canonical (no symbolic links, . or ..) path first
			p = canonical(p);
			// Append slash
			p /= "/";
			parametersFileDirectory = p.string(); 

			// Run error notification only after obtaining
			// all options and dealing with "--help"
			notify(vm);

			populateMaps(vm);

			// Set hardware default queue as required through provided options
			acl::hardware.setDefaultQueue(vm["platform"].as<string>(),
			                              vm["device"].as<string>());

			// Place it after(!) notify(vm);
			if (vm.count("check"))
			{
				cout << "Parameters consistency check - successful." << endl;
				exit(EXIT_SUCCESS);
			}
		}
		catch(exception& e)
		{
			errorMessage(string("ApplicationParametersManager::load() - ") + e.what());
		}
	}
		 
} //namespace asl