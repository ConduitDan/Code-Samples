#include "observers.h"

AreaObserver::AreaObserver(std::string path):FrameObserver(path) {
    file.open(_projectFolderPath +"area.csv",std::fstream::out);
}
AreaObserver::~AreaObserver() {
    file.close();
}
void AreaObserver::observe(Frame &f) {
    file<<f.area()<<',';
}
// VolumeObserver::VolumeObserver(std::string path):FrameObserver(path) {
//     file.open(_projectFolderPath +"volume.csv", std::fstream::out);
// }
// VolumeObserver::~VolumeObserver() {
//     file.close();
// }
// void VolumeObserver::observe(Frame &f) {
//     file<<f.volume()<<',';
// }

RadiusObserver::RadiusObserver(std::string path):FrameObserver(path) {
    file.open(_projectFolderPath +"radius.csv", std::fstream::out);
}
RadiusObserver::~RadiusObserver() {
    file.close();
}
void RadiusObserver::observe(Frame &f) {
    file<<f.radius()<<',';
}

SwapsObserver::SwapsObserver(std::string path):FrameObserver(path) {
    file.open(_projectFolderPath +"swaps.csv", std::fstream::out);
}
SwapsObserver::~SwapsObserver() {
    file.close();
}
void SwapsObserver::observe(Frame &f) {
    if (_pastFrame){
        file<<Geometry::bond_diff(_pastFrame->get_bonds(),f.get_bonds())<<',';
    }
    _pastFrame = &f;
}

VelocityProfileObserver::VelocityProfileObserver(std::string path):FrameObserver(path) {
    _nBins = ceil((_maxY-_minY)/_dy);
    _velocity = std::vector<double>(_nBins,0.0);
    _density = std::vector<unsigned int>(_nBins,0);
}

VelocityProfileObserver::VelocityProfileObserver(std::string path,int skip):FrameObserver(path) {
    _start = skip;
    _nBins = ceil((_maxY-_minY)/_dy);
    _velocity = std::vector<double>(_nBins,0.0);
    _density = std::vector<unsigned int>(_nBins,0);
}

void VelocityProfileObserver::write(){
    std::fstream file;
    file.open("velocity_profile.csv",std::fstream::out);
    double vel_ave = 0;
    for (int i = 0; i<_velocity.size(); i++){
        if (_density[i] == 0) vel_ave = 0;
        else vel_ave = _velocity[i]/_density[i];

        file<<_minY+_dy*i<<','<<vel_ave<<','<<_density[i]<<std::endl;
    }
    file.close();

}

void VelocityProfileObserver::observe(Frame &f) {
    // don't like this this math should be here not there on the frame. Frame is for holding data, but the computation here wehere it belongs
    if (_count>_start){
        f.vel_profile(_velocity,_density,_minY,_maxY,_dy);
    }
    _count++;
}
void MovieObserver::observe(Frame &f){
    std::string fileName = _projectFolderPath + "/images/" + std::to_string(_frame_counter);
    _povMaker->write_frame(&f,fileName);
    _frame_counter++;
}
MovieObserver::MovieObserver(std::string path):FrameObserver(path){
    _povMaker=std::make_unique<POVRayDrawer>();
}
ViscMovieObserver::ViscMovieObserver(std::string path, bool test):MovieObserver(path) {
    
    double left;
    double right;
    double top;
    double bottom;


    if (test) {
        left = -1;
        right = 5;
        top = 5;
        bottom =-1;
    } else {    
        left = 0;
        right = 100;
        top = 100;
        bottom = 0;
    }

    double centerx = (left + right)/2;
    double centery = (top + bottom)/2;
    
    double distance = 500;

    double largest = fmax(top-bottom,right-left)/2; 
    double angle = atan2(largest,distance)*2*360/(2*M_PI);
    _povMaker->set_viewangle(angle);
    _povMaker->set_viewpoint({centerx,centery,distance});
    _povMaker->set_look_at({centerx,centery,0});
}
void ViscMovieObserver::observe(Frame &f){
    MovieObserver::observe(f);
}

ForceObserver::ForceObserver(std::string force_type, std::string path, double k, double lmin, double lmax, int r):FrameObserver(path) {
    _params =std::make_unique<power_params>(power_params({lmax, lmin, k, r}));
    file.open(path + "/force_check.txt",std::fstream::out);
    if (force_type == "power")_force = Potentials::power;
    if (force_type == "power_harmonic")_force = Potentials::power_harmonic;
    
}

ForceObserver::~ForceObserver(){
    std::cout<<"pass is " <<_pass<<std::endl;
    file.close();
}
void ForceObserver::finish(){

    if (_pass){
        file<<"========================"<<std::endl;
        file<<"SUCCESS"<<std::endl;
        file<<"========================"<<std::endl;
    }
    file<<"Largest r = "<<_max_l<<std::endl;
    file<<"Smallest r = "<<_min_l<<std::endl;
}
void ForceObserver::observe(Frame &f){
    bool success = true;
    for (int nodeNum = 0; nodeNum<f.get_atoms().size(); nodeNum++){

    LammpsAtom node = f.get_atoms()[nodeNum];
    if (node.get_type()==2) continue;

    double fx = 0;
    double fy = 0;
    double fz = 0;
    
    std::array<double,3> dr;

    double dx, dy, dz;
    
    int other_node;
    std::unordered_set<int> bond_list=node.get_bonds();
    // gives us numbers of bonds that contain this atom
    for (auto &bond_no: bond_list){ 
        // find if the bond is starts from or ends at the node we want
        std::array<int,3> bond = f.get_bonds()[bond_no];
        if (bond[0] == nodeNum) { other_node = 1;}
        else if (bond[1] == nodeNum) { other_node = 0;}
        else throw "this bond doesn't contain the node we're looking for";

        dr = Geometry::minus(node.get_position(),f.get_atoms()[bond[other_node]].get_position());
        double r = sqrt(Geometry::dot(dr,dr));
        _max_l = fmax(_max_l,r);
        _min_l = fmin(_min_l,r);
        double f_scaled = _force(*_params,r);
        fx+=dr[0]*f_scaled;
        fy+=dr[1]*f_scaled;
        fz+=dr[2]*f_scaled;
    }
    std::array<double,3> diff = Geometry::minus(node.get_force(),std::array<double,3>({fx,fy,fz}));
    double error = sqrt(Geometry::dot(diff,diff));
    double fmag = sqrt(Geometry::dot(node.get_force(),node.get_force()));
    if (error/fmag>1e-2){
        file<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<std::endl;
        
    file<<"Checking Node "<<nodeNum<<" on frame "<< f.get_frame_number()<<" which clamis a force:"<<std::endl;
    file<<"fx: "<<node.get_force()[0]<<std::endl;
    file<<"fy: "<<node.get_force()[1]<<std::endl;
    file<<"fz: "<<node.get_force()[2]<<std::endl;
    
        file<<"I found that the force should be:"<<std::endl;
        file<<"fx: "<<fx<<std::endl;
        file<<"fy: "<<fy<<std::endl;
        file<<"fz: "<<fz<<std::endl;
        file<<"which is a error (relative) of: "<<error<<" ("<<error/fmag<<")"<<std::endl;
        for (auto &bond_no: bond_list){ 
            std::array<int,3> bond = f.get_bonds()[bond_no];
            if (bond[0] == nodeNum) { other_node = 1;}
            else if (bond[1] == nodeNum) { other_node = 0;}
            else throw "this bond doesn't contain the node we're looking for";

            dr = Geometry::minus(node.get_position(),f.get_atoms()[bond[other_node]].get_position());
            double r = sqrt(Geometry::dot(dr,dr));

            file<<"This node is connected to " <<bond[other_node]<<" with l="<<r<<" and dr = <"<<dr[0]<<','<<dr[1]<<','<<dr[2]<<">"<<std::endl;
            }

        success = false;
    }
    _pass *= success;  
}

    



}
