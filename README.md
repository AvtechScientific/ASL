
For more information, please visit <http://asl.org.il>.


# ASL

__Advanced Simulation Library (ASL)__ is a free and open source multiphysics simulation software package. Its computational engine is based, among others, on the [Lattice Boltzmann Methods](http://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) and is written in [OpenCL](http://en.wikipedia.org/wiki/OpenCL) which enable [extraordinarily efficient deployment](http://asl.org.il/benchmarks) on a variety of massively parallel architectures, ranging from inexpensive FPGAs, DSPs and GPUs up to heterogeneous clusters and supercomputers. The engine is hidden entirely behind C++ classes, so that no OpenCL knowledge is required from application programmers. ASL can be utilized to model various coupled physical and chemical phenomena and employed in a multitude of fields: computational fluid dynamics, virtual sensing, industrial process data validation and reconciliation, image-guided surgery, computer-aided engineering, high-performance scientific computing, crystallography, etc..


## License

ASL is distributed under the free GNU Affero General Public License (AGPLv3) with an optional [commercial license](http://asl.org.il/licensing).


## Support

Professional consulting, training and integration services are provided by [Avtech Scientific](http://avtechscientific.com), whose team created and continues to extend the library. The company offers [innovative R&D solutions](http://avtechscientific.com/services) and is involved in diverse academic and industrial [collaborative projects](http://avtechscientific.com/projects) dealing with complex multidisciplinary problems.


## Quick Start

### Installation

1. Install [cmake](http://cmake.org) (>=3.0.2, BSD License) and the required libraries:
	- [OpenCL](https://www.khronos.org/opencl) (>=1.2, OpenCL Specification License)
	- [boost](http://www.boost.org) (>=1.55, Boost Software License)
	- [VTK](http://vtk.org) (>=6.1, BSD License)
	- [optional](https://github.com/AvtechScientific/ASL/blob/master/cmake/ASLBuildOptions.cmake#L3): Matlab support with [matio](https://sourceforge.net/projects/matio) (>=1.5.2, BSD License)
	- [optional](https://github.com/AvtechScientific/ASL/blob/master/cmake/ASLBuildOptions.cmake#L4): API documentation with [doxygen](http://doxygen.org) (preferably with [graphviz](http://www.graphviz.org))
2. Download and extract the [ASL source code archive](https://github.com/AvtechScientific/ASL/releases/latest).
3. Create a build directory: `mkdir build-asl; cd build-asl`
4. Use [cmake generator](http://www.cmake.org/cmake/help/v3.2/manual/cmake-generators.7.html) to produce Makefiles: `cmake -G "Unix Makefiles" ../ASL` or project files for your IDE (Visual Studio, Xcode, Eclipse, etc.): `cmake -G "Visual Studio 10" ../ASL`
5. Run make (as root if installing into default destination `/usr/local`): `make install`

### Running an example

1. Go to examples: `cd examples/flow/locomotive`
2. Download geometry file [locomotive.stl](http://asl.org.il/input_data/locomotive.stl) from the [ASL input data page](http://asl.org.il/input_data).
3. Run: `./asl-locomotive --input locomotive.stl`  
Optionally: change parameters `./asl-locomotive --input locomotive.stl --dx 1 --dt 2` or write all of them into a file for later editing/reuse - `./asl-locomotive -g bigGrid.ini`. List all available options - `./asl-locomotive -h`.
4. Post-processing: see [step by step example](https://github.com/AvtechScientific/ASL/wiki/User-Guide#post-processing) and [locomotive.pvsm](http://asl.org.il/input_data/locomotive.pvsm) - the ParaView state file.

### Writing your own code using ASL

1. Take a look on [examples](http://asl.org.il/doc/Developer-Guide/examples.html) and the [API documentation](http://asl.org.il/doc/Developer-Guide/), start with [examples/flow/locomotive.cc](http://asl.org.il/doc/Developer-Guide/locomotive_8cc-example.html)
2. ASL installation supplies `ASL.pc` and `ASLConfig.cmake` files. To build your program using:

- `pkg-config`: ``c++ `pkg-config --cflags --libs ASL` -std=c++11 -o flow flow.cc``
- `cmake`: write a basic `CMakeLists.txt` file:

```cmake
project(locomotive)
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
find_package(ASL 0.1.4 CONFIG REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_executable(locomotive locomotive.cc)
target_link_libraries(locomotive PRIVATE ASL::aslnum ASL::aslvtk ASL::asl)
```