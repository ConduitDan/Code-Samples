#pragma once
#include <string>
#include <vector>
#include <stdio.h>
#include <array>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <exception>
#include <iostream>
#include <fstream>
#include <memory>
#include <limits>
#include "povInterface.h"
#include "forces.h"
#include "observers.h"

class POVRayDrawer;
struct hash_on_set_of_ints
: private std::hash<int>
{
  typedef int count_type;
  typedef std::hash<count_type> base;
  std::size_t operator()(std::unordered_set<int> const&obj) const
  {
    return base::operator()(std::accumulate(obj.begin(),obj.end(),count_type()));
  }
};



class FrameObserver;
class Frame;
class LammpsAtom;
class LammpsPostProcessor {
protected: 
    std::string _projectFolderPath;
    std::vector<Frame> _trajectoryData;
public:
    LammpsPostProcessor(std::string folder): _projectFolderPath(folder){};
    void read_data(int start_ts, int final_ts, int increment);
    void calculate_area();
    void calculate_radius();
    void calculate_swaps();
    void make_movie();
    void calcV();
    void make_movie_visc();
    void write_frame_to_input();
    bool check_flips();
    void post_process_frames(std::tuple<int,int,int> frameInfo, std::vector<FrameObserver*> frameObs);
};

class Frame {
protected:
    std::vector<LammpsAtom> _atoms;
    std::vector<std::array<int,3>> _bonds; // bond is {start, end, type}
    std::vector<std::array<int,3>> _facets;
    std::vector<std::array<std::array<int,4>,5>> _flips_start;
    std::vector<std::array<std::array<int,4>,5>> _flips_end;
    std::vector<int> _connectivity;

    int _frameNumber;
public:
    Frame(){};
    Frame(int frame_no): _frameNumber(frame_no){}; 
    bool readAtom(const char* filename);
    bool readFacet(const char* filename);
    bool readBond(const char* filename);
    bool read_and_check_Flips(const char* filename,std::string projectFolder);
    bool check_dihedral_substitution(int pre[4], int post[4], int subs[2]);
    bool check_dihedral_connectivity(int dihedral[4]);
    bool check_valid_dihedral(int dihedral[4]);
    void generate_connectivity();

    void vel_profile(std::vector<double> &velocity, std::vector<unsigned int> &density,double minY,double maxY, double dy);

    double area();
    double radius();
    std::vector<LammpsAtom>& get_atoms() {return _atoms;}
    std::vector<std::array<int,3>>& get_bonds(){return _bonds;}
    int get_frame_number(){return _frameNumber;}
};


class LammpsAtom {
protected:
    std::array<double,3> _position;
    std::array<double,3> _velocity;
    std::array<double,3> _force;
    std::unordered_set<int> _bonds;
    int _type;
    int _mol;
public:
    LammpsAtom(){}
    LammpsAtom(int type_in,int mol_in,std::array<double,3> r,std::array<double,3> v, std::array<double,3> f):_type(type_in),_mol(mol_in),_position(r),_velocity(v),_force(f){}
    std::array<double,3> get_position(){return _position;}
    std::array<double,3> get_velocity(){return _velocity;}
    std::array<double,3> get_force(){return _force;}

    int get_type(){return _type;}
    bool add_bond(int bond){auto in = _bonds.insert(bond); return in.second;}
    bool remove_bond(int bond){return _bonds.erase(bond);}
    std::unordered_set<int> get_bonds(){return _bonds;}

};

