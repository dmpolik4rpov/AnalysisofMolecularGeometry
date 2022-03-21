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
	cout << "Число атомов: " << mol.natom << endl;
	cout << "Координаты атомов в прямоугольной системе координат:\n";
	mol.print_geom();



	cout << "Межатомные расстояния (Ангстремы):\n";
	for (int i = 0; i < mol.natom; i++)
			for(int j = 0; j<i; j++)
				printf("%d %d %8.5f\n", i, j, mol.bond(i, j));



	cout << "\nУглы связей:\n";
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