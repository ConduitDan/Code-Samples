#pragma once
#include "forces.h"
#include "lammpsPostProcessor.h"
#include "povInterface.h"
class POVRayDrawer;
class Frame;


class FrameObserver {
protected:
std::string _projectFolderPath;

public:
FrameObserver(std::string path):_projectFolderPath{path}{}

virtual void observe(Frame &f) = 0;
virtual void finish(){}
};

class AreaObserver: public FrameObserver {
protected:
    std::fstream file;
public:
    AreaObserver(std::string path);
    ~AreaObserver();
    void observe(Frame &f);
};
// class VolumeObserver: public FrameObserver {
// protected:
//     std::fstream file;
// public:
//     VolumeObserver(std::string path);
//     ~VolumeObserver();
//     void observe(Frame &f);
// };
class VelocityProfileObserver: public FrameObserver {
protected:
    double _dy = 0.1;
    double _minY = 25;
    double _maxY = 75;
    int _nBins;
    std::vector<double> _velocity;
    std::vector<unsigned int> _density;
    int _count = 0;
    int _start = 0;

public:
    VelocityProfileObserver(std::string path);
    VelocityProfileObserver(std::string path,int skip);
    void write();
    void observe(Frame &f);
    void finish(){std::cout<<"writing to file in: "<<_projectFolderPath<<std::endl; write();}
};

class MovieObserver: public FrameObserver {
protected:
    std::unique_ptr<POVRayDrawer> _povMaker;
    long _frame_counter = 0;
public:
    MovieObserver(std::string path);

    virtual void observe(Frame &f);
};
class ViscMovieObserver: public MovieObserver{
public:
    ViscMovieObserver(std::string path,bool test);
    void observe(Frame &f);
};

class SwapsObserver: public FrameObserver{
protected:
    Frame *_pastFrame = nullptr;
    std::fstream file;

public:
    SwapsObserver(std::string path);
    ~SwapsObserver();
    void observe(Frame &f);

};
class RadiusObserver: public FrameObserver{
protected:
    std::fstream file;
public:

    RadiusObserver(std::string path);
    ~RadiusObserver();
    void observe(Frame &f);

};

struct power_params;
class ForceObserver: public FrameObserver{
protected:
    std::unique_ptr<power_params> _params;
    std::fstream file;
    double (*_force)(power_params,double);
    int _nodeNum;
    bool _pass = true;
    double _max_l = 0;
    double _min_l = 1e4;
public:
    ForceObserver(std::string force_type, std::string path, double k, double lmin, double lmax, int r);
    ~ForceObserver();
    void finish();
    void observe(Frame &f);
};
