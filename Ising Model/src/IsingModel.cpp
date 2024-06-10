#include "IsingModel.h"

void MetropolisHasting::run(int MCSteps, IsingModel* model) {
    double deltaE;
    IsingSpin *testSpin;
    for (int i = 0; i<MCSteps; i++){
        for (int j = 0; j <model->get_system_size(); j++){
            // claculate the energy
            testSpin = &(model->get_random_spin());
            // here we use that E(sigma) = -E(-sigma); I.E our hamiltion for E_i is linear in sigma_i
            // delta E is E_i(-sigma_i)-E_I(sigma_i) = -2E_i(sigma_i)
            deltaE = -2*model->get_energy(*testSpin);
            if (accept(deltaE)) {
                testSpin->flip();
            }

        }
    }
}

bool MetropolisHasting::accept(double deltaE) { 

    // if this reduces energy accept the change
    if (deltaE<0) {
        return true; 
    } else {
        // otherwise generate a random number and compare exp(-delta /T to it)
        double r = _distr(_gen);
        if (r<exp(-deltaE/(_T+1e-16))) {
            return true;
        }
    }
    return false;
}

IsingSpin& IsingModel::get_random_spin() {
    int i = distr(_gen);
    return (*_grid)[i];
}

double IsingModel::get_energy(IsingSpin& spin) {
    double E = 0;
    int index = spin.get_index();

    for (_neighbors->begin(index); !_neighbors->is_done(); _neighbors->next()){
        E-= _J*spin.get_state() * (**_neighbors).get_state();
    }
    
    E-= spin.get_state()*_H;
    return E;
}

IsingModel::IsingModel(int L, double H, double J,std::string latticeType): _L{L},_H{H},_J{J} {
    if (latticeType == "Square") {
        _grid = std::make_unique<SquareGrid<IsingSpin>>(L,L);

    } else if (latticeType =="Triangular") {
        _grid = std::make_unique<TriangularGrid<IsingSpin>>(L,L);
        
    } else {
        throw std::domain_error("Grid type must be either \"Square\" or \"Triangular\"");
    }

    _neighbors = std::unique_ptr<GridNeighborIterator<IsingSpin>>(_grid->get_neighbor_iterator());
    
    for (int i = 0; i< _L*_L; i++){
        (*_grid)[i].set_index(i);
    }
    distr=std::uniform_int_distribution<>(0,_L*_L);
}

void IsingExperiment::phase_diagram(int L,std::string latticeType) {
    IsingModel model(L,0,1,latticeType);
    MetropolisHasting MCAlgorithm(10.0);
    MagnetismIsingObserver obs("phaseDiagram.out");
    // 101 points from 10 to 0
    for (double T = 10; T>=0; T-=0.1){
        MCAlgorithm.set_T(T);
        MCAlgorithm.run(100,&model);
        obs.takeMeasurement(&model,&MCAlgorithm,T);
    }

}

void IsingExperiment::hysteresis(int L,std::string latticeType) {
    IsingModel model(L,0,1,latticeType);
    MetropolisHasting MCAlgorithm(10);
    MagnetismIsingObserver obs("hysteresisT=10.out");
    // 101 points from 0 to 10
    for (double H = -10; H<=10; H+=0.1){
        model.set_H(H);
        MCAlgorithm.run(100,&model);
        obs.takeMeasurement(&model,&MCAlgorithm,H);
    }
    for (double H = 10; H>=-10; H-=0.1){
        model.set_H(H);
        MCAlgorithm.run(100,&model);
        obs.takeMeasurement(&model,&MCAlgorithm,H);
    }

    MagnetismIsingObserver obs2("hysteresisT=0.5.out");
    MCAlgorithm.set_T(0.5);
    for (double H = -10; H<=10; H+=0.1){
        model.set_H(H);
        MCAlgorithm.run(100,&model);
        obs2.takeMeasurement(&model,&MCAlgorithm,H);
    }
    for (double H = 10; H>=-10; H-=0.1){
        model.set_H(H);
        MCAlgorithm.run(100,&model);
        obs2.takeMeasurement(&model,&MCAlgorithm,H);
    }
}
void MagnetismIsingObserver::takeMeasurement(IsingModel *subject, MetropolisHasting *MCAlgorithm ,double value){
    double m = 0;
    for (int j = 0; j<_numMeasurements; j++){
        for (int i = 0; i<subject->get_system_size(); i++){
            m+=subject->get_spin(i).get_state();
        }
        MCAlgorithm->run(1,subject);
    }
    m=m/subject->get_system_size()/_numMeasurements;
    outFile<<value<<","<<m<<std::endl;
}

void MovieIsingObserver::takeMeasurement(IsingModel *subject, MetropolisHasting *MCAlgorithm ,double value){
    outFile.open("data/"+_fileName+std::to_string(_counter)+".out",std::fstream::out);
    outFile<<value<<',';
    std::array<double,2> coord; 
    for (int i = 0; i<subject->get_system_size(); i++){
        coord = subject->get_coord(i);
        outFile<<coord[0]<<','<<coord[1]<<','<<subject->get_spin(i).get_state()<<',';
    }
    outFile.close();
    _counter++;
}

void IsingExperiment::anneal_movie(int L,std::string latticeType) {
    IsingModel model(L,0,1,latticeType);
    MetropolisHasting MCAlgorithm(10.0);
    MovieIsingObserver obs("IsingMovie"+latticeType+".out");
    // 101 points from 10 to 0
    for (double T = 5; T>=0; T-=0.1){
        MCAlgorithm.set_T(T);
        for (int i = 0; i<30;i++){
            
            MCAlgorithm.run(1,&model);
            obs.takeMeasurement(&model,&MCAlgorithm,T);

        }
        std::cout<<(5-T)*20<<'%'<<std::endl;
    }
}
