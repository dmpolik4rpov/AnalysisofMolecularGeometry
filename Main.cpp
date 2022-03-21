#include "Molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;

int main()
{
	setlocale(LC_ALL, "Russian");

	Molecule mol("geom.dat", 0);
	cout << "Number of atoms: " << mol.natom << endl;
	cout << "Atom coordinates in Cartesian coord system:\n";
	mol.print_geom();



	cout << "Interatomic distances in Angstrem:\n";
	for (int i = 0; i < mol.natom; i++)
			for(int j = 0; j<i; j++)
				printf("%d %d %8.5f\n", i, j, mol.bond(i, j));



	cout << "\Bond angle:\n";
	for (int i = 0; i < mol.natom; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < j; k++) {
				if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0)
					printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i, j, k) * (180.0 / acos(-1.0)));
			}
		}
	}





	return 0;
}
