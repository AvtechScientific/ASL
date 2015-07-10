# ASL

__Advanced Simulation Library (ASL)__ is a free and open source multiphysics simulation software package. Its computational engine is based, among others, on the [Lattice Boltzmann Methods](http://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) and is written in [OpenCL](http://en.wikipedia.org/wiki/OpenCL) which enable [extraordinarily efficient deployment](http://asl.org.il/benchmarks) on a variety of massively parallel architectures, ranging from inexpensive FPGAs, DSPs and GPUs up to heterogeneous clusters and supercomputers. The engine is hidden entirely behind C++ classes, so that no OpenCL knowledge is required from application programmers. ASL can be utilized to model various coupled physical and chemical phenomena and employed in a multitude of fields: computational fluid dynamics, virtual sensing, industrial process data validation and reconciliation, image-guided surgery, computer-aided engineering, high-performance scientific computing, etc..


## License

ASL is distributed under the free GNU Affero General Public License (AGPLv3) with an optional [commercial license](http://asl.org.il/licensing). Professional support and consulting services are provided by [Avtech Scientific](http://avtechscientific.com), whose team created and continues to extend the library. The company offers [innovative R&D solutions and services](http://avtechscientific.com/services) and is involved in diverse academic and industrial [collaborative projects](http://avtechscientific.com/projects) dealing with complex multidisciplinary problems.


## Further information

For more information, please visit <http://asl.org.il>.


## Quick Start

### Installation

1. Install [cmake](http://cmake.org) (BSD License) and the required libraries:
	- [OpenCL](https://www.khronos.org/opencl) (OpenCL Specification License)
	- [C++ bindings for OpenCL](https://www.khronos.org/registry/cl/api/1.1/cl.hpp) (OpenCL Specification License)
	- [boost](http://www.boost.org) (Boost Software License)
	- [VTK](http://vtk.org) (BSD License)
	- optional: [matio](https://sourceforge.net/projects/matio) (BSD License)
	- optional: [doxygen](http://doxygen.org) (preferably with [graphviz](http://www.graphviz.org)) for API documentation generation
2. Download and extract the [ASL source code archive](https://github.com/AvtechScientific/ASL/releases/latest).
3. Create a build directory: `mkdir build-asl; cd build-asl`
4. Use [cmake generator](http://www.cmake.org/cmake/help/v3.2/manual/cmake-generators.7.html) to produce Makefiles: `cmake -G "Unix Makefiles" ../ASL` or project files for your IDE (Visual Studio, Xcode, Eclipse, etc.): `cmake -G "Visual Studio 10" ../ASL`
5. Run make (as root if installing into default destination `/usr/local`): `make install`

### Running an example

1. Go to examples: `cd examples/flow/locomotive_in_tunnel`
2. Copy the .stl input file: `cp ../../../../ASL/examples/input_data/locomotive.stl .`
3. Run: `./locomotive_in_tunnel`. Optionally: change some parameters - `./locomotive_in_tunnel --dx 0.1 --dt 2` or write all of them into a file for later editing/reuse - `./locomotive_in_tunnel -g bigGrid.ini`. List all available options - `locomotive_in_tunnel -h`.
4. Post-processing: see [step by step example](https://github.com/AvtechScientific/ASL/wiki/User-Guide#post-processing) and the state file `examples/input_data/locomotive_in_tunnel.pvsm`.

### Writing your own code using ASL

1. Take a look on examples, start with [examples/flow/locomotive_in_tunnel.cc](http://asl.org.il/doc/Developer-Guide/locomotive_in_tunnel_8cc-example.html)
2. ASL installation supplies `ASL.pc` and `ASLConfig.cmake` files. To build your program using:
	- `pkg-config` - launch ``c++ `pkg-config --cflags --libs ASL` -o locomotive_in_tunnel locomotive_in_tunnel.cc``
	- `cmake` - write a basic `CMakeLists.txt` file:

```
project(locomotive)
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
find_package(ASL 0.1.4 CONFIG REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_executable(locomotive_in_tunnel locomotive_in_tunnel.cc)
target_link_libraries(locomotive_in_tunnel PRIVATE ASL::aslnum ASL::aslvtk ASL::asl)
```
