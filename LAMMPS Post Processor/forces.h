#pragma once
#include "lammpsPostProcessor.h"

class LammpsAtom;

class Geometry {

public:
    static double area(std::vector<LammpsAtom> vertices, std::vector<std::array<int,3>> facets);
    static double radius(std::vector<LammpsAtom> vertices,int type);
    static std::array<double,3> cross(std::array<double,3> v1, std::array<double,3> v2);
    static std::array<double,3> minus(std::array<double,3> v1, std::array<double,3> v2);
    static double dot(std::array<double,3> v1, std::array<double,3> v2);
    static int bond_diff(std::vector<std::array<int,3>> bonds1,std::vector<std::array<int,3>> bonds2);
    static double length(LammpsAtom v1, LammpsAtom v2);


};
struct sw_harmonic_params{
    double rmin;
    double rmax;
    double delta;
    double sigma;
    double leq;
    double k;
};
struct power_params{
    double lmax;
    double lmin;
    double k;
    int r;
};

class Potentials {
public:

    static double sw_harmonic(sw_harmonic_params params, std::array<int,2> bond,std::vector<LammpsAtom>* atoms);
    static double power_harmonic(power_params params, double l);
    static double power(power_params params, double l);
    static double dihedral(double stiffness, std::array<int,4> dihedral,std::vector<LammpsAtom>* atoms);
};

