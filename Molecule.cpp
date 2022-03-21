#include "Molecule.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>



Molecule mol("geom.dat", 0);

void Molecule::print_geom()
{
    for (int i = 0; i < natom; i++)
        printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

void Molecule::translate(double x, double y, double z)
{
    for (int i = 0; i < natom; i++) {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}


Molecule::Molecule(const char* filename, int q)
{
    charge = q;

    // Readin the file
    std::ifstream is(filename);
    assert(is.good());

    // Counting the number of atoms
    is >> natom;

    // New memory
    zvals = new int[natom];
    geom = new double* [natom];
    for (int i = 0; i < natom; i++)
        geom[i] = new double[3];

    for (int i = 0; i < natom; i++)
        is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

    is.close();
}

Molecule::~Molecule()
{
    delete[] zvals;
    for (int i = 0; i < natom; i++)
        delete[] geom[i];
    delete[] geom;
}


double Molecule::bond(int a, int b) {

    double dangstr = sqrt((geom[a][0] - geom[b][0]) * (geom[a][0] - geom[b][0]) + (geom[a][1] - geom[b][1]) * (geom[a][1] - geom[b][1]) + (geom[a][2] - geom[b][2]) * (geom[a][2] - geom[b][2]));
        dangstr = dangstr * 1.88973;///now in angstrem
        return dangstr;
}
double Molecule::unit(int cart, int a, int b)
{
    return -(geom[a][cart] - geom[b][cart]) / bond(a, b);
}
double Molecule::angle(int a, int b, int c)
{
    return acos(unit(0, b, a) * unit(0, b, c) + unit(1, b, a) * unit(1, b, c) + unit(2, b, a) * unit(2, b, c));
}

