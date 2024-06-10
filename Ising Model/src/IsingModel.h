#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <iostream>
#include <memory>
#include "Grid.h"



// class SquareGridIterator{
// protected:
//     int _L;
//     IsingModel *parent;
// public:
//     SquareGridNeighborIterator get_neighbor_iterator(int i);
//     int get_L(){return _L;}

// };


class IsingSpin;

class IsingModel{
protected:
    // SquareGrid _grid;
    double _J;
    double _H;
    int _L;
    std::mt19937 _gen;
    std::uniform_int_distribution<> distr;

    std::unique_ptr<Grid<IsingSpin>> _grid;
    std::unique_ptr<GridNeighborIterator<IsingSpin>> _neighbors;


    int boundaryCheck(int index);


public:
    IsingModel(int L,double H,double J,std::string lattice_type);

    IsingSpin& get_random_spin();
    IsingSpin& get_spin(int i) {return (*_grid)[i];};
    
    double get_energy(IsingSpin& spin);
    

    void set_seed(double seed) {_gen = std::mt19937(seed);}
    void set_H(double H){_H = H;}
    int get_system_size(){return _grid->get_num_elements();}

    std::array<double,2> get_coord(int i){return _grid->get_coord(i);}

};

class IsingSpin {
protected:
    bool _state = false;
    int _index;
public:
    int get_state(){if (_state) return 1; else return -1;}
    int get_index(){return _index;}

    void flip(){_state = !_state;}
    void set_index(int i){_index = i;}
};

class MetropolisHasting {
protected:
    double _T;
    std::mt19937 _gen;
    std::uniform_real_distribution<> _distr;
    
public:
    MetropolisHasting(double T): _T{T}{ set_seed(1); _distr=std::uniform_real_distribution<>(0,1);}
    void run (int MCSteps,IsingModel* model);
    void set_T(double T){_T=T;}
    bool accept(double deltaE);
    void set_seed(double seed) {_gen = std::mt19937(seed);}

};
class IsingObserver{
protected:
    std::fstream outFile;
    int _numMeasurements = 100;
public:
    IsingObserver(){}
    IsingObserver(std::string fileName){
        outFile.open(fileName,std::fstream::out);
    }
    ~IsingObserver(){
        outFile.close();
    }
    virtual void takeMeasurement(IsingModel *subject,MetropolisHasting* MC, double value){

    }

};

class MagnetismIsingObserver: public IsingObserver{
public:
    using IsingObserver::IsingObserver;
    void takeMeasurement(IsingModel *subject, MetropolisHasting* MC, double value);
};
class MovieIsingObserver: public IsingObserver{
protected:
    std::string _fileName;
    int _counter=0;
public:
    MovieIsingObserver(std::string fileName){
        _fileName=_fileName;
    }
    void takeMeasurement(IsingModel *subject, MetropolisHasting* MC, double value);
};



class IsingExperiment {
public:
    static void phase_diagram(int L, std::string latticeType);
    static void hysteresis(int L, std::string latticeType);
    static void anneal_movie(int L, std::string latticeType);

};