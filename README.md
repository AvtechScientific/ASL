# ASL

__Advanced Simulation Library (ASL)__ is a free and open source multiphysics simulation software package. Its computational engine is based, among others, on the [Lattice Boltzmann Methods](http://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) and is written in [OpenCL](http://en.wikipedia.org/wiki/OpenCL) which enable [extraordinarily efficient deployment](http://asl.org.il/benchmarks) on a variety of massively parallel architectures, ranging from inexpensive FPGAs, DSPs and GPUs up to heterogeneous clusters and supercomputers. The engine is hidden entirely behind C++ classes, so that no OpenCL knowledge is required from application programmers. ASL can be utilized to model various coupled physical and chemical phenomena and employed in a multitude of fields: computational fluid dynamics, virtual sensing, industrial process data validation and reconciliation, image-guided surgery, computer-aided engineering, high-performance scientific computing, etc..


## License

ASL is distributed under the free GNU Affero General Public License (AGPLv3) with an optional [commercial license](http://asl.org.il/licensing). Professional support and consulting services are provided by [Avtech Scientific](http://avtechscientific.com), whose team created and continues to extend the library. The company offers [innovative R&D solutions and services](http://avtechscientific.com/services) and is involved in diverse academic and industrial [collaborative projects](http://avtechscientific.com/projects) dealing with complex multidisciplinary problems.


## Quick Start

### Installation

1. Install [cmake](http://cmake.org) (BSD License) and the required libraries:
	- [OpenCL](https://www.khronos.org/opencl) (OpenCL Specification License)
	- [C++ bindings for OpenCL](https://www.khronos.org/registry/cl/api/1.1/cl.hpp) (OpenCL Specification License)
	- [boost](http://www.boost.org) (Boost Software License)
	- [VTK](http://vtk.org) (BSD License)
	- optional: [matio](https://sourceforge.net/projects/matio) (BSD License)
2. Download and extract the [ASL source code archive](https://github.com/AvtechScientific/ASL/releases/latest).
3. Create a build directory: `mkdir build-asl; cd build-asl`
4. Use [cmake generator](http://www.cmake.org/cmake/help/v3.2/manual/cmake-generators.7.html) to produce Makefiles: `cmake -G "Unix Makefiles" ../ASL` or project files for your IDE (Visual Studio, Xcode, Eclipse, etc.): `cmake -G "Visual Studio 10" ../ASL`
5. Run make (as root if installing into default destination `/usr/local`): `make install`

### Running an example

1. Go to tests: `cd examples/flow/locomotive_in_tunnel`
2. Copy the .stl input file: `cp ../../../../ASL/examples/input_data/locomotive.stl .`
3. Run: `./locomotive_in_tunnel`. Optionally: change some parameters - `./locomotive_in_tunnel --dx 0.1 --dt 2` or write all of them into a file for later editing/reuse - `./locomotive_in_tunnel -g bigGrid.ini`. See `locomotive_in_tunnel -h` for more information.
4. Post-processing: [step by step example](https://github.com/AvtechScientific/ASL/wiki/User-Guide#post-processing).

### Writing your own code using ASL

1. Take a look on examples, e.g. `examples/flow/locomotive_in_tunnel.cc`
2. To build your program using `cmake` see e.g. `examples/flow/CMakeLists.txt`
3. To build your program with tools others than `cmake` run `make VERBOSE=1` and/or consult `CMakeCache.txt` to get better understanding of the compiler flags, library paths and dependencies involved. The output of `make install` shows the location of installed public include headers and libraries (by default: `/usr/local/include/asl-X.Y.Z` and `/usr/local/lib`).


## Further information

For more information, please visit <http://asl.org.il>.
