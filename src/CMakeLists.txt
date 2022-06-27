cmake_minimum_required(VERSION 3.2)	 
			
project(AnalisysofMolecularGeometry)

set(CMAKE_BUILD_TYPE "Release")	

set(SOURCE_EXE Main.cpp Molecule.cpp)					

add_executable(Main ${SOURCE_EXE})	

target_link_libraries(Main Molecule)		