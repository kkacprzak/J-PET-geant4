# Installation

## Required packages:  
* `cmake`  
* `root 6.0`  
* [`geant4.10.04`](https://github.com/Geant4/geant4)  
 Built with cmake flags:  
 `-DGEANT4_USE_QT=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_GDML=ON`  
* [cadmesh](https://github.com/christopherpoole/CADMesh)
* doxygen

## How to compile?
`mkdir build`  
`cd build`  
`cmake .. -DCMAKE_PREFIX_PATH=path_to_the_directory_containing_cadmesh-config.cmake`  
`make`  
output file: (in build folder)  
`bin/jpet_mc`  

## How to create documentation?
(in build folder)  
`cmake .. && make doc`  
- open the `doc/html/index.html` in your favorite web browser  

