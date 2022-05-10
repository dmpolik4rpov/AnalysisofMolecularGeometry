#include "Molecule.h"


int main()
{

	Molecule mol("geom.dat", 0);
	std::cout << "Number of atoms:" << mol.natom << std::endl;
	std::cout << "Input Cartesian coordinates : \n";
	mol.printGeom();



	std::cout << "Interatomic distances (Bohr):\n";
	for (int i = 0; i < mol.natom; i++)
		for (int j = 0; j < i; j++)
			std::printf("%d %d %8.5f\n", i, j, mol.bond(i, j));



	std::cout << "\nBond angles (degrees):\n";
	for (int i = 0; i < mol.natom; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++) {
				if (mol.bond(i, j) < 4 && mol.bond(j, k) < 4)
					std::printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i, j, k) * (180.0 / std::acos(-1.0)));
			}
		}
	}


	std::cout << "\nOut of Plane (degrees):\n";

	for (int i = 0; i < mol.natom; i++) {
		for (int k = 0; k < mol.natom; k++) {
			for (int j = 0; j < mol.natom; j++) {
				for (int l = 0; l < j; l++) {
					if (i != j && i != k && i != l && j != k && k != l && mol.bond(i, k) < 4.0 && mol.bond(k, j) < 4.0 && mol.bond(k, l) < 4.0)
						std::printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.outOfPlane(i, j, k, l) * (180.0 / acos(-1.0)));
				}
			}
		}
	}
	mol.PrintCOM();
		return 0;
}
