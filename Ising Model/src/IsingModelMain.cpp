#include "IsingModel.h"


typedef enum{PhaseDiagram,Hysteresis,AnnealMovie} ExperimentType;
typedef enum{Square, Triangular} LatticeType;

int main(int argc, char** argv){

    int entry;


    std::cout<<"Choose the experiment"<<std::endl;
    std::cout<<"Enter 0 for Phase Diagram, 1 for Hysteresis and 2 for Annealing Movie"<<std::endl;
    std::cin>>entry;
    ExperimentType experiment = static_cast<ExperimentType>(entry);

    int systemSize;
    std::cout<<"Chose a system size"<<std::endl;
    std::cin>>systemSize;

    // parameter validation
    if (systemSize<0) throw std::invalid_argument( "System Size must be positive" );

    std::cout<<"Choose a lattice type"<<std::endl;
    std::cout<<"Enter 0 for Square, 1 for Triangular"<<std::endl;
    std::cin>>entry;

    LatticeType myLattice = static_cast<LatticeType>(entry);

    std::string latticeType;
    if (myLattice == Square) latticeType = "Square";
    if (myLattice == Triangular) latticeType ="Triangular";
    


    // choose the experiment and run it
    if (experiment == PhaseDiagram) {
        IsingExperiment::phase_diagram(systemSize,latticeType);
    } else if (experiment == Hysteresis) {
        IsingExperiment::hysteresis(systemSize,latticeType);
    } else if (experiment == AnnealMovie) {
        IsingExperiment::anneal_movie(systemSize,latticeType);
    } 

    return 0;
}