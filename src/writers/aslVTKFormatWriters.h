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


#ifndef ASLVTKFORMATWRITERS_H
#define ASLVTKFORMATWRITERS_H

#include "aslWriter.h"
//#include <data/aslDataWrapper.h>

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

namespace asl 
{
	/// Writes data accumulated in Writer into a file with VTK XML format "vti"
	/// \ingroup IO
	class WriterVTKXML : public Writer
	{
		Block *newBl;
		public:
			WriterVTKXML(const std::string & file, Block *nbl=0);
			void write();
	};


	/// writes \p data in a file with VTK XML format "vti"
	/// \ingroup IO
	void writeVTKXML(const std::string & fileName,
	                 const AbstractData & data,
	                 const std::string & name);

} // asl

#endif // ASLVTKFORMATWRITERS_H

