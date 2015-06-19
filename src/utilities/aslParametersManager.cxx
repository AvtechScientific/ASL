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
#include <math/aslVectorsDynamicLength.h>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <acl/aclTypesList.h>
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
	Parameter<T>::Parameter(string key_,
	                        string description_,
	                        string units_):
		key(key_),
		description(description_),
		units(units_)
	{
		// Adds itself to current ParametersManager
		if (ParametersManager::current == NULL)
			errorMessage("ParametersManager was not instantiated and is not available");

		WildcardCheck<T>::check(key, parameter.p);

		ParametersManager::current->add(parameter,
		                                key,
		                                description);
	}

	template Parameter<string>::Parameter(string key,
	                                      string description,
	                                      string units);
		
	#define BOOST_TT_rep_expression(r, data, t) \
		template Parameter<t>::Parameter(string key, \
										 string description, \
										 string units); \
		template Parameter<AVec<t>>::Parameter(string key, \
											   string description, \
											   string units); \
		template Parameter<map<string, t>>::Parameter(string key, \
													  string description, \
													  string units); \
		template Parameter<map<string, AVec<t>>>::Parameter(string key, \
															string description, \
															string units);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression


	template <typename T> Parameter<T>::Parameter(T defaultValue,
	                                              string key_,
	                                              string description_,
	                                              string units_):
		key(key_),
		description(description_),
		units(units_)
	{
		// Adds itself to current ParametersManager
		if (ParametersManager::current == NULL)
			errorMessage("ParametersManager was not instantiated and is not available");
		
		ParametersManager::current->add(parameter,
		                                defaultValue,
		                                key,
		                                description);
	}

	template Parameter<string>::Parameter(string defaultValue,
	                                      string key,
	                                      string description,
	                                      string units);
		
	#define BOOST_TT_rep_expression(r, data, t) \
		template Parameter<t>::Parameter(t defaultValue, \
										 string key, \
										 string description, \
										 string units); \
		template Parameter<AVec<t>>::Parameter(AVec<t> defaultValue, \
											   string key, \
											   string description, \
											   string units);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression


	ParametersManager * ParametersManager::current(NULL);

	ParametersManager::ParametersManager():
		configurationOptions("Configuration options"),
		parametersFileStr(""),
		folder(""),
		folderWithSlash(""),
		programName(""),
		programVersion("")
	{
		enable();
		// Add platform and device parameters with default values
		Parameter<string> platform(acl::getPlatformVendor(acl::hardware.defaultQueue),
									"platform", "Default computation platform", "");
		Parameter<string> device(acl::getDeviceName(acl::hardware.defaultQueue),
									"device", "Default computation device", "");
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
	                                                  string key,
	                                                  string description)
	{
		configurationOptions.add_options()
			(key.c_str(), value<T>(&parameter.v())->required(), description.c_str());

		// Add the option to the default parameters file			
		parametersFileStr += "\n# " + description + "\n" + key + " = \n";
	}


	template<typename T> void ParametersManager::add(UValue<map<string, T>> parameter,
	                                                  string key,
	                                                  string description)
	{
		configurationOptions.add_options()
			(key.c_str(), value<T>()->required(), description.c_str());
		// ToDo: maps - Add the option to the default parameters file
//		parametersFileStr += "\n# " + description + "\n" + key + " = " + numToStr(defaultValue) + "\n";
	}


	template<typename T> void ParametersManager::add(UValue<T> parameter,
	                                                 T defaultValue,
	                                                 string key,
	                                                 string description)
	{		
		configurationOptions.add_options()
			(key.c_str(), value<T>(&parameter.v())->default_value(defaultValue), description.c_str());

		// Add the option to the default parameters file
		parametersFileStr += "\n# " + description + "\n" + key + " = " + numToStr(defaultValue) + "\n";
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


	void ParametersManager::load(int argc, char * argv[],
	                             string programName_,
	                             string programVersion_)
	{
		programName = programName_;
		programVersion = programVersion_;

		variables_map vm;

		options_description genericOptions("Generic options");

		genericOptions.add_options()
			("help,h", "display this help and exit")
			("version,v", "display version and exit")
			("devices,d", "display available devices and exit")
			("folder,f", value<string>()->default_value("Default"),
			 "path to the working folder that contains configuration file - parameters.ini")
			("parameters,p", value<string>()->default_value("."),
			 "generate default configuration file parameters.ini, write it to the provided path and exit")
			("check,c", "check configuration for consistency and exit");

		positional_options_description positional;
		positional.add("folder", 1);

		options_description allOptions;

		allOptions.add(genericOptions).add(configurationOptions);

		try
		{
			store(command_line_parser(argc, argv).options(allOptions).positional(positional).run(), vm);

			if (vm.count("help"))
			{
				cout << "Usage: " + programName + " [WORKING_FOLDER] [OPTION]...\n"
					 << allOptions
					 << endl;
				exit(0);
			}

			if (vm.count("version"))
			{
				cout << programName + " " + programVersion
					 << endl;
				exit(0);
			}

			if (vm.count("devices"))
			{
				cout << programName + " " + programVersion + "\n\n"
					<< "Default computation device:"
					<< acl::hardware.getDefaultDeviceInfo() << "\n\n"
					<< "List of all available platforms and their devices:\n"
					<< acl::hardware.getDevicesInfo()
					<< endl;
				exit(0);
			}

			if (vm.count("parameters"))
			{
				path p(vm["parameters"].as<string>());
				// add at least one slash at the end
				p /= "/";
				// and then cut all possible slashes at the end
				p = p.parent_path();
				p /= "/";
				p /= "parameters.ini";

				cout << "Writing default configuration file to: "
					 << p.string() << endl;

				writeParametersFile(p.string());
				exit(0);
			}

			path p(vm["folder"].as<string>());
			// add at least one slash at the end
			p /= "/";
			// and then cut all possible slashes at the end
			p = p.parent_path();
			folder = p.string();
			p /= "/";
			folderWithSlash = p.string();
			p /= "parameters.ini";
			ifstream ifs(p.string());
			if (!ifs)
			{
				// Only warn, since all options might have default values, or required values
				// provided through the command line - so no configuration file is required
				warningMessage("ParametersManager::load() - can not open configuration file: " + p.string());
			}

			parsed_options parsed = parse_config_file(ifs, allOptions, true);
			store(parsed, vm);
			// Run error notification only after obtaining
			// all options and dealing with "--help"
			notify(vm);

			populateMaps(vm);

			// Set Hardware default queue as required through provided options
			acl::hardware.setDefaultQueue(vm["platform"].as<string>(), vm["device"].as<string>());

			// Place it after(!) notify(vm);
			if (vm.count("check"))
			{
				cout << programName + " " + programVersion + "\n"
					 << "Consistency check - successful."
					 << endl;
				exit(0);
			}
		}
		catch(exception& e)
		{
			errorMessage(string("ParametersManager::load() - ") + e.what());
		}
	}


	void ParametersManager::load(string configFile)
	{
		variables_map vm;

		try
		{
			ifstream ifs(configFile);
			if (!ifs)
				errorMessage("Can not open configuration file: " + configFile);

			parsed_options parsed = parse_config_file(ifs, configurationOptions, true);
			store(parsed, vm);
			notify(vm);
			populateMaps(vm);
		}
		catch(exception& e)
		{
			errorMessage(string("ParametersManager::load() - ") + e.what());
		}
	}


	string ParametersManager::getFolder()
	{
		return folder;
	}


	string ParametersManager::getFolderWithSlash()
	{
		return folderWithSlash;
	}
	

	void ParametersManager::writeParametersFile(const std::string fileName)
	{
		ofstream fo(fileName);
		if (!fo)
			errorMessage("ParametersManager::writeParametersFile() - can not open file: " + fileName);
			
		// Prepend informative header
		parametersFileStr = "# Parameters file with default values (where available).\n# Generated by "
									+ programName + " version " + programVersion + "\n\n"
									+ "# Get the list of all available computation devices by running:\n"
									+ "# `" + programName + " -d`\n" + parametersFileStr;
									
		fo << parametersFileStr;
		fo.close();
	}

} //namespace asl
