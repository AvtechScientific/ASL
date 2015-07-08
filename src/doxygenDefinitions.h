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


/**
	\mainpage
	\dotfile mainDiagram.dot 
*/


// This file contains Doxygen definitions


/// \defgroup ACL Advanced Computational Language


/// \defgroup Operators Operators
/// \ingroup ACL


/// \defgroup Functions Functions
/// \ingroup ACL


/// \defgroup MathOperators Mathematical Operators
/// \ingroup Operators


/**
 \defgroup AssignmentOperators Assignment Operators
 \ingroup Operators
*/

/**
 \defgroup BooleanOperators Boolean Operators
 \ingroup Operators
*/

/**
 \defgroup MathFunctions Mathematical Functions
 \ingroup Functions
*/

/**
 \defgroup SpecialPurposeFunctions Special Purpose Functions
 \ingroup Functions
*/

/**
 \defgroup SynchronizationFunctions Synchronization Functions
 \ingroup Functions
*/


/// \defgroup HostInterectionFunctions Host Interaction Functions
/// \ingroup Functions


/**
 \defgroup ControlStructures Control Structures
 \ingroup ACL
*/

/**
 \defgroup KernelGen OpenCL Kernel Generation
 \ingroup ACL
*/

/**
 \defgroup ComplexDataTypes Complex Data Types
 \ingroup ACL
*/


/// \defgroup HardwareInformation Hardware Information
/// \ingroup ACL

/// \defgroup AtomicFunctions Atomic Functions
/// \ingroup ACL


/**
\defgroup generateVE VectorOfElements generators
 \related VectorOfElements
 \ingroup ComplexDataTypes
*/

/**
\defgroup generateME MatrixOfElements generators
 \related MatrixOfElements
 \ingroup ComplexDataTypes
*/


/**
 \defgroup Numerics Numerics
 */

/**
 \defgroup NumMethods Numerical Methods
 \ingroup Numerics
 */

/**
 \defgroup BoundaryConditions Boundary Conditions
 \ingroup Numerics
*/


/**
 \defgroup Physics Physics
*/


/**
 \defgroup TransportProcesses Transport Processes
 \ingroup Physics
*/

/**
 \defgroup ChemicalReactions Chemical Reactions
 \ingroup Physics
*/

/**
 \defgroup GenericBC Generic Boundary Conditions
 \ingroup BoundaryConditions
*/


/**
 \defgroup TransportProcessesBC Transport Processes Boundary Conditions
 \ingroup TransportProcesses
 \ingroup BoundaryConditions
 */

/**
 \defgroup ElasticityBC Elasticity Boundary Conditions
 \ingroup Elasticity
 \ingroup BoundaryConditions
 */


/// \defgroup Interfacing Interfacing with other software


/**
 \defgroup IO Input/Output
 \ingroup Interfacing
 */


/// \defgroup DataAnalysis Data Analysis Tools
/// \ingroup Numerics


/// \defgroup DataUtilities Data utilities
/// \ingroup Numerics

/// \defgroup LevelSet Advanced Level Set Method
/// \ingroup Numerics

/// \defgroup PF Position Function
/// \ingroup Numerics


/// \defgroup Geom Geometric primitives and their manipulations


/// \defgroup LDI Library Design Issues


/// \defgroup Utilities Utilities


/// \defgroup ErrorMessaging Error Messaging
/// \ingroup Utilities


/**
  \defgroup Splines Splines
*/



namespace boost
{ 
	template<class T> class shared_ptr {}; 
	template<class T> class shared_array {}; 
	template<class T> class scoped_ptr {}; 
	template<class T> class scoped_array {}; 
}

namespace std
{ 
	template<class T> class vector {}; 
	template<class T> class shared_ptr {}; 
	template<class T> class shared_array {}; 
	template<class T> class unique_ptr {}; 
	template<class T> class unique_array {}; 
} 
