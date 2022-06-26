#include "Molecule.h"



Molecule mol("geom.dat", 0);

void Molecule::printGeom()
{
    for (int i = 0; i < natom; i++)
       std::printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
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
    std::ifstream is(filename);
    is >> natom;
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
    return dangstr;
}
///in the axes direction(cart = 0 = x, cart = 1 = y, cart = 2 = z)
double Molecule::unit(int cart, int a, int b)
{
    return -(geom[a][cart] - geom[b][cart]) / bond(a, b);
}
double Molecule::angle(int a, int b, int c)
{
   double AngRad = std::acos(unit(0, b, a) * unit(0, b, c) + unit(1, b, a) * unit(1, b, c) + unit(2, b, a) * unit(2, b, c));
   return AngRad;///here calculation comes in Radians however in Main.cpp there is a conversion into degrees
}

double Molecule::outOfPlane (int a, int b, int c, int d) {


 double ebcd_x = (mol.unit(1, c, b) * mol.unit(2, c, d) - mol.unit(2, c, b) * mol.unit(1, c, d));
 double ebcd_y = (mol.unit(2, c, b) * mol.unit(0, c, d) - mol.unit(0, c, b) * mol.unit(2, c, d));
 double ebcd_z = (mol.unit(0, c, b) * mol.unit(1, c, d) - mol.unit(1, c, b) * mol.unit(0, c, d));


  double ex_x = ebcd_x * mol.unit(0, c, a);
  double ey_y = ebcd_y * mol.unit(1, c, a);
  double ez_z = ebcd_z * mol.unit(2, c, a);

  double phi = mol.angle(b,c,d);

  double theta = (ex_x + ey_y + ez_z) / std::sin(phi);
  if (theta < -1.0)
  {
      theta = std::asin(-1.0);
  }
  else if (theta > 1.0) {
      theta = std::asin(1.0);
  }
  else theta = std::asin(theta);

    return theta;
}
void Molecule::PrintCOM() {
    double Masses[135];
    std::ifstream input("masses.txt");
    for (int i = 0; i < 33; i++) {////at this moment masses.txt contains only 33 elements. Going to add more in the future
        input >> Masses[i];
    };
    double M = 0.0;
    for (int i = 0; i < mol.natom; i++) {
        M += Masses[mol.zvals[i]];///calculating the molar mass
    }
    double xcm = 0.0;///
    double ycm = 0.0;
    double zcm = 0.0;
    double mi;
    for (int i = 0; i < mol.natom; i++) {
        mi = Masses[mol.zvals[i]];
        xcm += mi * mol.geom[i][0];///calculating the centre of masses on every asis
        ycm += mi * mol.geom[i][1];
        zcm += mi * mol.geom[i][2];
    }
    xcm /= M;
    ycm /= M;
    zcm /= M;
    printf("\nMolecular center of mass: %12.8f %12.8f %12.8f\n", xcm, ycm, zcm);

    mol.translate(-xcm, -ycm, -zcm);

    std::cout << "Molar mass is " << M << ":\n";
}
