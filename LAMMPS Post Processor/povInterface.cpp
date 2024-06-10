#include "povInterface.h"

void POVRayDrawer::write_frame(Frame *f, std::string imageName) {
    std::string povFile= imageName+".pov";
    file.open(povFile, std::fstream::out );
    if (!file.is_open()){
        fprintf(stderr, "Could not open pov file %s \n",povFile.c_str());
        return;
    }
    write_header();
    // for each atom in the vertex write it in the pov file
    

    
    for (auto &bond :f->get_bonds()){
        write_bond(f->get_atoms()[bond[0]],f->get_atoms()[bond[1]]);
    }

    int i = 0;
    for (auto &vertex :f->get_atoms()){
        write_vertex(vertex);
	//write_text(std::to_string(i),vertex.get_position());
	i++;
    }
    file.close();
    std::string povCommand = "povray ./\""+imageName+".pov\" -D +A +W"+_width+" +H"+_height+" All_Console=Off > /dev/null";
    system(povCommand.c_str());

}

std::string POVRayDrawer::pov_vector(double *v) {
    std::string out = "<" + std::to_string(v[0])+", "\
                     + std::to_string(v[1])+", "\
                     + std::to_string(v[2])+">";
    return out;
}

// std::string POVRayDrawer::color() { return std::string("<0.5, 0.5, 0.5>"); }

void POVRayDrawer::write_header() {

    file<<"#include \"colors.inc\"";
    file<<"background { rgb "+_background+ "}";
    file<<"camera {"\
        <<"location "<<_viewpoint\
        <<"up <0,1,0> right <-1.33,0,0> angle "<<_viewangle\
        <<"look_at " << _look_at <<" sky " << _sky << " }"<<std::endl;

    for (auto &light: _lights) {
        file<<"light_source {"<<light <<" color White}"<<std::endl;
    }
}

void POVRayDrawer::write_vertex(LammpsAtom vertex) {
    std::string out = "sphere { "
              + pov_vector(vertex.get_position().data())+ ", "
              + std::to_string(_vertex_radius) 
              + " texture { pigment { rgb " + _vertex_color
              + "} } }";
    file<<out<<std::endl;
    return;
}

void POVRayDrawer::write_bond(LammpsAtom v1, LammpsAtom v2) {
    double length = Geometry::length(v1,v2);
    // for viscosity, ignore giant bonds
    if (length > 50) {return;}
    

    double stress = abs(length-1);
    double color[3];
    double stressmax = 1.67;
    color[0] = stress/stressmax;
    color[1] = 0;
    color[2] = (stressmax-stress)/stressmax;
    
    color[0] = color[0]/sqrt(color[0]*color[0]+color[2]*color[2]);
    color[2] = color[2]/sqrt(color[0]*color[0]+color[2]*color[2]);


    double radius = 0.5*length*_bondAspectratio;
    std::string out = "cylinder { "\
                     + pov_vector(v1.get_position().data()) + ", "\
                     + pov_vector(v2.get_position().data()) + ", "\
                     + std::to_string(radius)\
                    + " texture { pigment { rgb "+ pov_vector(color) + " } } }";
    file<<out<<std::endl;
    return;
}
void POVRayDrawer::write_text(std::string txt, std::array<double,3> placement){
std::string out = "text { ttf \"timrom.ttf\"";

	out+= " \""+txt + "\"0.1, 0 pigment{White} " ;
	out+= "scale <0.35, 0.35, 0.35> ";
	out+= "translate " + std::to_string(placement[0]) + "*x " ;
	out+= "translate " + std::to_string(placement[1]) + "*y}" ;
	file<<out<<std::endl;
}
