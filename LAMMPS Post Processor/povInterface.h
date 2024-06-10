#ifndef POV_INTERFACE_H
#define POV_INTERFACE_H
#include "lammpsPostProcessor.h"
//#include "observers.h"
#include "forces.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
class Frame;
class LammpsAtom;


class POVRayDrawer{
private:
    std::fstream file;
    double _bondAspectratio = 0.1;
    double _vertex_radius = 0.15;
    
    std::string _bond_color = "<0.5, 0.5, 0.5>";
    std::string _vertex_color = "<0.2, 0.75, 0.2>";
    std::vector<std::string> _lights = {"<40, 40, 150>"};//, "<-50, 50, 200>"};
    std::string _viewangle = "24";
    std::string _look_at = "<0, 0, 0>";
    std::string _sky = "<0, 1, 0>";
    std::string _viewpoint = "<0, 0, 200>";
    std::string _width = "2048";//"4096";
    std::string _height = "1536";//"";3072
    std::string _background = "<0.05, 0.05, 0.05>";//"<1, 1, 1>"

    /*    self.camera = camera
    if (!self.camera) self.camera = Camera() // Default Camera
    self.graphic = graphic

    self.antialias = self.camera.antialias
    self.width = self.camera.width
    self.height = self.camera.height
    self.viewangle = self.camera.viewangle
    self.viewpoint = self.camera.viewpoint
    self.look_at = self.camera.look_at
    self.sky = self.camera.sky
    self.light = nil*/

    std::string pov_vector(double* v);
    



public:

    void set_viewangle(double angle){_viewangle = std::to_string(angle);}
    void set_look_at(std::array<double,3> point){_look_at = pov_vector(point.data());}
    void set_viewpoint(std::array<double,3> point){_viewpoint = pov_vector(point.data());}
    void set_background(std::array<double,3> point){_background = pov_vector(point.data());}

    POVRayDrawer(){}
    void write_frame(Frame *f,std::string imageName);
    void write_header();
    void write_vertex(LammpsAtom vertex);
    void write_bond(LammpsAtom v1, LammpsAtom v2);
    void write_text(std::string txt, std::array<double,3> placement);
    void render();
};


#endif
