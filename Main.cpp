#include "Molecule.h"


int main()
{

	Molecule mol("geom.dat", 0);
	std::cout << "Number of atoms:" << mol.natom << std::endl;
	std::cout << "Input Cartesian coordinates : \n";
	mol.print_geom();



	std::cout << "Interatomic distances (Bohr):\n";
	for (int i = 0; i < mol.natom; i++)
			for(int j = 0; j<i; j++)
				std::printf("%d %d %8.5f\n", i, j, mol.bond(i, j));



	std::cout << "\nBond angles:\n";
	///cout << mol.natom;
	for (int i = 0; i < mol.natom; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++) {
				if (mol.bond(i, j) < 4 && mol.bond(j, k) < 4)
					std::printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i, j, k) * (180.0 / std::acos(-1.0)));
			}
		}
	}





	return 0;
}
