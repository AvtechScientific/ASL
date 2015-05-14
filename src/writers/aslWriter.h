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


#ifndef ASLWRITER_H
#define ASLWRITER_H

//#include <acl/aclMath/aclVectorOfElementsDef.h>
#include<string>
#include<memory>
#include<vector>

namespace acl
{
	class VectorOfElements;
	class VectorOfElementsData;
	class Kernel;
}

namespace asl 
{
	class AbstractData;
	class Block;

	class Writer
	{
		public:
			Writer(const std::string & file_);
			~Writer();
			virtual void write() = 0;
			void enable();
			void addScalars(std::string name, AbstractData & data);
			void addVector(std::string name, AbstractData & data);
			void addScalars(std::string name, acl::VectorOfElementsData & data);
			void addVector(std::string name, acl::VectorOfElementsData & data);
			void addScalars(std::string name, 
			                const acl::VectorOfElements & data, 
			                acl::Kernel & kernel, 
			                unsigned int nGhost = 1);
			void addVector(std::string name,
			               const acl::VectorOfElements & data, 
			               acl::Kernel & kernel,
			               unsigned int nGhost = 1);

			static Writer * current;

		protected:
			std::shared_ptr<Block> block;
			std::vector<std::pair<std::string, acl::VectorOfElementsData>> scalarFields;
			std::vector<std::pair<std::string, acl::VectorOfElementsData>> vectorFields;
			std::string file;
			unsigned int numOfWrites;
	};

} // asl

#endif // ASLWRITER_H
