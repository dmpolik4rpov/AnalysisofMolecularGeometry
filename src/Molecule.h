#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>




class Molecule
{
public:
	int natom;///number of atoms
	int charge;
	int* zvals;///number of protons 
	double** geom;
	std::string point_group;

	void printGeom();
	void translate(double x, double y, double z);
	double bond(int atom1, int atom2);
	double angle(int atom1, int atom2, int atom3);
	double outOfPlane(int atom1, int atom2, int atom3, int atom4);
	double unit(int cart, int atom1, int atom2);
	void PrintCOM();

	Molecule(const char* filename, int q);
	~Molecule();
};


