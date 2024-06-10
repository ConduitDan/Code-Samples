#include "lammpsPostProcessor.h"


void LammpsPostProcessor::post_process_frames(std::tuple<int, int, int> frameInfo,
                                              std::vector<FrameObserver*> frameObs) {
    std::string atoms_file;
    std::string bond_file;
    std::string dihedral_file;
    Frame newFrame;
    int start_ts = std::get<0>(frameInfo);
    int final_ts = std::get<1>(frameInfo);
    int increment = std::get<2>(frameInfo);
    
    if (start_ts == final_ts) increment=1;

    bool success = true;
    //#pragma omp parallel for private(atoms_file,bond_file,dihedral_file,newFrame,success)
    for (int frameNumber = start_ts; frameNumber<=final_ts; frameNumber+=increment){
        atoms_file = _projectFolderPath + "/trajectory/"+ std::string(std::to_string(frameNumber)) + ".dump";
        bond_file = _projectFolderPath + "/topology/bond"+ std::string(std::to_string(frameNumber)) + ".dump";
        dihedral_file = _projectFolderPath + "/topology/facet"+ std::string(std::to_string(frameNumber)) + ".dump";
        newFrame = Frame(frameNumber);
        success *= newFrame.readAtom(atoms_file.c_str()); 
        success *= newFrame.readBond(bond_file.c_str()); 
        newFrame.readFacet(dihedral_file.c_str());
        if (!success) {
            std::cout<<"Failed to read dump "<<frameNumber<<std::endl;
            _trajectoryData.resize((increment-start_ts)/increment);
            return;
        }
        for (auto &obs: frameObs){
            obs->observe(newFrame);
        }
    }
}

bool Frame::readAtom(const char* fileName) {
    std::ifstream fp;

    // char line[200];
    std::string fileLine;

    fp.open(fileName);
    if (!fp.is_open()) {
        perror("Could not open atom file");
        std::cout<<"atom file: "<<fileName<<std::endl;
        // printf(stderr, fileName);
        // printf(stderr, "\n");
        
        return false;
    }
    // skip the first 3 lines
    for (int i = 0; i < 3; i++) std::getline(fp,fileLine);//fgets(line, 50, fp);  // read line n

    // read the number of atoms
    int numAtoms;
    std::getline(fp,fileLine);
    sscanf(fileLine.c_str(), "%d\n", &numAtoms);
    _atoms.resize(numAtoms);

    // skip the next 5 lines
    for (int i = 0; i < 5; i++) std::getline(fp,fileLine);//fgets(line, 50, fp);  // read line n

    // read in the id mol type x y z vx vy vz
    int id, type, mol;
    std::array<double,3> r;
    std::array<double,3> v;
    std::array<double,3> f = {0,0,0};
    for (int i = 0; i < numAtoms; i++) {
        std::getline(fp,fileLine);
        sscanf(fileLine.c_str(), "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &id, &type, &mol, &r[0], &r[1], &r[2], &v[0], &v[1], &v[2],&f[0],&f[1],&f[2]);
        //index by 0
        id--;

        if (id<numAtoms && id>=0) {
            _atoms[id] = LammpsAtom(type, mol, r, v,f);
        } else {
            std::cout<<"Tried using ID "<<id<< " but number of atoms is "<<numAtoms<<std::endl;
            std::cout<<"Line is '"<<fileLine<<"'"<<std::endl;
            throw "bad id number read atoms";
        }
    }
    fp.close();
    return true;
}

bool Frame::readBond(const char* fileName) {
    FILE* fp;
    char line[200];
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        perror("Could not open bond file");
        std::cout<<"bond file: "<<fileName<<std::endl;

        // fprintf(stderr, fileName);
        // fprintf(stderr, "\n");
        return false;
    }
    // skip the first 3 lines
    for (int i = 0; i < 3; i++) fgets(line, 50, fp);  // read line n
    // read the number of bonds
    int numBonds;
    fscanf(fp, "%d\n", &numBonds);
    _bonds.resize(numBonds);

    // skip 5
    for (int i = 0; i < 5; i++) fgets(line, 50, fp);  // read line n

    int id, type, source, target;
    for (int i = 0; i < numBonds; i++) {
        fscanf(fp,"%d %d %d %d\n", &id, &type, &source, &target);
        // index by 0
        source--;
        target--;

        _bonds[i] = {source, target, type};
        _atoms[source].add_bond(i);
        _atoms[target].add_bond(i);
    }
    
    fclose(fp);
    return true;
}
bool Frame::read_and_check_Flips(const char* filename,std::string projectFolderPath) {


    // ok so the order they come out corresponds to the map ID but the bonds correspond to the tag ID
    // I can work with this. I guess I can construct my own damn tag->map thing


    // current connundrum

    std::fstream fp;
    POVRayDrawer povMaker;
    //{rmin, rmax, delta, sigma, leq, k}

    sw_harmonic_params sw_params = {0.1,2.67,0.5,10000,1,100};
    double di_stiffness = 500;
    int high_energy_counter = 0;
    double highE;
    // FILE* fp;
    fp.open(filename,std::ios::in);
    std::string line;
    if (!fp.is_open()) {
        perror("Could not open flips file");
        // printf(stderr, fileName);
        // printf(stderr, "\n");
        return false;
    }

    int a,b,c,d,a1,b1,c1,d1;
    int neighbors_old[4][4];
    int neighbors_new[4][4];
    int neighbors_reps[4][2];

//main: {8719, 9086, 9087, 9104} -> {9086, 8719, 9104, 9087};  Neighbors:\
//       {9087, 9086, 8719, 8717} -> {9104, 9086, 8719, 8717} 9087->9104,\
//       {9086, 9087, 8719, 9080} -> {9104, 9087, 8719, 9080} 9086->9104,\
//       {9086, 9087, 9104, 9100} -> {8719, 9087, 9104, 9100} 9086->8719,\
//       {9087, 9086, 9104, 9102} -> {8719, 9086, 9104, 9102} 9087->8719
    std::string Parse = "main: {%d,%d,%d,%d} -> {%d,%d,%d,%d}; Neighbors: ";
    std::string neighbor_str = "{%d,%d,%d,%d} -> {%d,%d,%d,%d} %d->%d";
    Parse = Parse + neighbor_str + ',' + neighbor_str +',' + neighbor_str +',' + neighbor_str + " Energy: %lf" +'\n';
    std::string positionParse = "{<%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf>}\n";
    double atomPos[4][3];
    generate_connectivity();
    int count = 0;
    double E;
    while (getline(fp,line)){
        sscanf(line.c_str(), Parse.c_str(),\
         &a,&b,&c,&d,&a1,&b1,&c1,&d1,\
         &neighbors_old[0][0],&neighbors_old[0][1],&neighbors_old[0][2],&neighbors_old[0][3],\
         &neighbors_new[0][0],&neighbors_new[0][1],&neighbors_new[0][2],&neighbors_new[0][3],\
         &neighbors_reps[0][0],&neighbors_reps[0][1],\

         &neighbors_old[1][0],&neighbors_old[1][1],&neighbors_old[1][2],&neighbors_old[1][3],\
         &neighbors_new[1][0],&neighbors_new[1][1],&neighbors_new[1][2],&neighbors_new[1][3],\
         &neighbors_reps[1][0],&neighbors_reps[1][1],\

         &neighbors_old[2][0],&neighbors_old[2][1],&neighbors_old[2][2],&neighbors_old[2][3],\
         &neighbors_new[2][0],&neighbors_new[2][1],&neighbors_new[2][2],&neighbors_new[2][3],\
         &neighbors_reps[2][0],&neighbors_reps[2][1],\

         &neighbors_old[3][0],&neighbors_old[3][1],&neighbors_old[3][2],&neighbors_old[3][3],\
         &neighbors_new[3][0],&neighbors_new[3][1],&neighbors_new[3][2],&neighbors_new[3][3],\
         &neighbors_reps[3][0],&neighbors_reps[3][1],&E);

        getline(fp,line);
        sscanf(line.c_str(), positionParse.c_str(),\
            &atomPos[0][0],&atomPos[0][1],&atomPos[0][2],\
            &atomPos[1][0],&atomPos[1][1],&atomPos[1][2],\
            &atomPos[2][0],&atomPos[2][1],&atomPos[2][2],\
            &atomPos[3][0],&atomPos[3][1],&atomPos[3][2]);



        std::array<std::array<int,4>,5> finalDihedrals;
        finalDihedrals[0] = {a1,b1,c1,d1};
        finalDihedrals[1] = {neighbors_new[0][0],neighbors_new[0][1],neighbors_new[0][2],neighbors_new[0][3]};
        finalDihedrals[2] = {neighbors_new[1][0],neighbors_new[1][1],neighbors_new[1][2],neighbors_new[1][3]};
        finalDihedrals[3] = {neighbors_new[2][0],neighbors_new[2][1],neighbors_new[2][2],neighbors_new[2][3]};
        finalDihedrals[4] = {neighbors_new[3][0],neighbors_new[3][1],neighbors_new[3][2],neighbors_new[3][3]};

        std::array<std::array<int,4>,5> initalDihedrals;
        initalDihedrals[0] = {a,b,c,d};
        initalDihedrals[1] = {neighbors_old[0][0],neighbors_old[0][1],neighbors_old[0][2],neighbors_old[0][3]};
        initalDihedrals[2] = {neighbors_old[1][0],neighbors_old[1][1],neighbors_old[1][2],neighbors_old[1][3]};
        initalDihedrals[3] = {neighbors_old[2][0],neighbors_old[2][1],neighbors_old[2][2],neighbors_old[2][3]};
        initalDihedrals[4] = {neighbors_old[3][0],neighbors_old[3][1],neighbors_old[3][2],neighbors_old[3][3]};

        bool success = true;
        for (int neighbor_number = 1; neighbor_number<5; neighbor_number++) {
            // check substituion validity
            success = success && check_dihedral_substitution(initalDihedrals[neighbor_number].data(),finalDihedrals[neighbor_number].data(),neighbors_reps[neighbor_number-1]);

        }
        // check connectivity
        int dihedral[4]={a,b,c,d};
        success = success &&check_dihedral_connectivity(dihedral);
        //check bonds and dihedral corespondence
        success = success && check_valid_dihedral(dihedral);
        if (!success){
            std::cout<<"Failed on line: "<<count<<std::endl<<line<<std::endl;
            return false;
        }

        double delta=0;
        delta += Potentials::sw_harmonic(sw_params,{a,d},&_atoms);
        if (!std::isinf(delta)) delta -= Potentials::sw_harmonic(sw_params,{b,c},&_atoms);

        for (int neighbor_number = 0; neighbor_number<5; neighbor_number++) {
            if (!std::isinf(delta)) delta+= Potentials::dihedral(di_stiffness,finalDihedrals[neighbor_number],&_atoms);
            if (!std::isinf(delta)) delta-= Potentials::dihedral(di_stiffness,initalDihedrals[neighbor_number],&_atoms);

        }

        if (delta>0) {
            if (exp(-delta)<1e-4) {
                high_energy_counter++;
                highE = delta;
                // std::cout<<"Found High Energy flip "<<delta<<" In {"\
                // <<a<<", "<<b<<", "<<c<<", "<<d<<"} -> {"\
                // <<b<<", "<<a<<", "<<d<<", "<<c<<"}"<<std::endl;


            }
        }
        if (std::isinf(delta)) return false;
        if (abs((delta-E)/std::max(delta,E))>1) {
            std::cout<<" Wrong energy calculated. Lammps Reports "<<E<< " by my caluclation it should be "<<delta<<std::endl;
            std::cout<<" In {"\
                <<a<<", "<<b<<", "<<c<<", "<<d<<"} -> {"\
                <<b<<", "<<a<<", "<<d<<", "<<c<<"}"<<std::endl;
            return false;
        }


        std::string fileName = projectFolderPath + "/images/" + std::to_string(count);
        // povMaker.write_frame(this,fileName);



        count++;
    }
    fp.close();
    if ((double)high_energy_counter/(double)count > 1e-2) {
        std::cout<<"Too many High energy jumps, "<< high_energy_counter<<" out of "<< count\
        <<" had a probability of 1e-4 or less"<<std::endl;
        std::cout<<"The last high energy jump had an energy of "<<highE<<std::endl;
        return false;
    } 


    return true;

}
bool Frame::check_dihedral_substitution(int pre[4], int post[4], int subs[2]) { 
    // make sure pre has the from and not the to
    // make sure post has the to and not the from
    // make sure the rest match

    bool passed;
    bool found_to = false;
    bool found_from = false;
    for (int i = 0; i<4; i++) {
        if (pre[i] == subs[0]) found_from = true;
        if (pre[i] == subs[1]) found_to = true;
    }
    // test passes if we have found the from but not the to
    passed = found_from && !found_to;

    if (!passed) {
        std::cout<<"Bad neighbor dihedral, {"<<pre[0]<<", "<<pre[1]<<", "<<pre[2]<<", "<<pre[3]<<"}"<<std::endl;
        std::cout<<"Wants to make substitution "<< subs[0] << " -> "<<subs[1]<<std::endl;
        return false;
    }
    found_to = false;
    found_from = false;
    for (int i = 0; i<4; i++) {
        if (post[i] == subs[0]) found_from = true;
        if (post[i] == subs[1]) found_to = true;
    }
    // test passes if we have found the to but not the from
    passed = !found_from && found_to;
    if (!passed) {
        std::cout<<"Bad neighbor dihedral, {"<<post[0]<<", "<<post[1]<<", "<<post[2]<<", "<<post[3]<<"}"<<std::endl;
        std::cout<<"Tried to make substitution "<< subs[0] << " -> "<<subs[1]<<std::endl;
        return false;
    }
    return true;

}
bool Frame::check_dihedral_connectivity(int dihedral[4]) { 
    
    //we have to subrtact 1 from all the node numbers for indexing issues

    // for both b and c make sure connectiviy is between 4 and 9 and at the start and finish
    auto test_bc = [&](int node) {
        if(_connectivity[node]<4) {
            std::cout<<"Deficient node "<<node+1<<" with connectivity "<< _connectivity[node]<<" in swap"<<std::endl;
            return false;
        }
        return true;
    };
    auto test_ad = [&](int node) {
        if (_connectivity[node]>9) {
                std::cout<<"Overfull node "<<node+1<<" with connectivity "<< _connectivity[node]<<" in swap"<<std::endl;
                std::cout<<"Bonds on this node are:"<<std::endl;
                for (auto &bond :_atoms[node].get_bonds()){
                    std::cout<<"{"<<_bonds[bond][0]+1<<","<<_bonds[bond][1]+1<<"}"<<std::endl;


                }
                std::cout<<"Plus the new one for this flip"<<std::endl;
            return false;

        }
        return true;
    };
    bool a = test_ad(dihedral[0]-1);
    bool b = test_bc(dihedral[1]-1);
    bool c = test_bc(dihedral[2]-1);
    bool d = test_ad(dihedral[3]-1);
    bool success = (a && b && c && d);
    if (!success) return false;
    _connectivity[dihedral[0]-1]++;
    _connectivity[dihedral[1]-1]--;
    _connectivity[dihedral[2]-1]--;
    _connectivity[dihedral[3]-1]++;
    a = test_ad(dihedral[0]-1);
    b = test_bc(dihedral[1]-1);
    c = test_bc(dihedral[2]-1);
    d = test_ad(dihedral[3]-1);
    success = (a && b && c && d);
    if (!success) return false;
    return true;

}
bool Frame::check_valid_dihedral(int dihedral[4]) { 
//  alpha ---- a ---- beta    alpha ---- a ---- beta
//       \    / \    /             \    /|\    /
//        \  1   2  /               \  1 | 2  /
//         \/     \/                 \/  |  \/
//          b -0- c       ===>        b  0  c
//         /\     /\                 /\  |  /\       
//        /  3   4  \               /  3 | 4  \      
//       /    \ /    \             /    \|/    \     
//  gamma ---- d ---- delta   gamma ---- d ---- delta
// if this dihedral is to be valid then we need to find the bonds
// a-b, b-d, d-c, c-a, b-c;
// then we should swap b-c->d-a;
    int bonds_to_find[5][2];
    bonds_to_find[0][0] = dihedral[0];bonds_to_find[0][1] = dihedral[1]; //a-b
    bonds_to_find[1][0] = dihedral[1];bonds_to_find[1][1] = dihedral[3]; //b-d
    bonds_to_find[2][0] = dihedral[3];bonds_to_find[2][1] = dihedral[2]; //d-c
    bonds_to_find[3][0] = dihedral[2];bonds_to_find[3][1] = dihedral[0]; //c-a
    bonds_to_find[4][0] = dihedral[1];bonds_to_find[4][1] = dihedral[2]; //b-c
        
    // function that finds the connection from node1 to node2
    auto find_bond = [&](int node1, int node2){
        // remeber that bonds are off by one because lamps indexs by 1 and we index by 0
        node1--;
        node2--;
        std::unordered_set<int> node1bonds = _atoms[node1].get_bonds();
        // std::cout<<"checking bonds on node "<<node1<<std::endl;
        for (auto &bond_num:node1bonds){
            if (_bonds[bond_num][0]==node2 ||_bonds[bond_num][1]==node2){
                return bond_num;
            }
            // std::cout<<"{"<<_bonds[bond_num][0]<<","<<_bonds[bond_num][1]<<"}"<<std::endl;
        }
        return -1;
    };

    int bond_num1;
    int bond_num2;

    for (int i = 0; i<5; i++) {
        bond_num1 = find_bond(bonds_to_find[i][0],bonds_to_find[i][1]);
        bond_num2 = find_bond(bonds_to_find[i][1],bonds_to_find[i][0]);

        if (bond_num1!=bond_num2 || bond_num1==-1) { 
            std::cout<<"dihedral {"<<dihedral[0]<<','<<dihedral[1]<<','<<dihedral[2]<<','<<dihedral[3]\
            <<"} does not have supporting bonds, offending bond is {"\
            <<bonds_to_find[i][0]<<','<<bonds_to_find[i][1]<<"} bond nums are: "\
            <<bond_num1<<" and "<<bond_num2<<std::endl;
            return false;
        }

    }
    // now flip b-c to a-d

    _bonds[bond_num1][0]=dihedral[0]-1;
    _bonds[bond_num1][1]=dihedral[3]-1;
    
    _atoms[dihedral[1]-1].remove_bond(bond_num1);
    _atoms[dihedral[2]-1].remove_bond(bond_num1);

    _atoms[dihedral[0]-1].add_bond(bond_num1);
    _atoms[dihedral[3]-1].add_bond(bond_num1);
    return true;
}
void Frame::generate_connectivity() {
    // connectivity is the size of atoms
    _connectivity.resize(_atoms.size());
    // fill it with 0s
    std::fill(_connectivity.begin(), _connectivity.end(), 0);

    // loop over bonds
    for (int count = 0; count<_bonds.size(); count++ ) {
        auto bond = _bonds[count];
        // if a node partisipates in a bond it has +1 connectivity
        _connectivity[bond[0]]++;
        _connectivity[bond[1]]++;
        _atoms[bond[0]].add_bond(count);
        _atoms[bond[1]].add_bond(count);
    }
}
void Frame::vel_profile(std::vector<double>& velocity, std::vector<unsigned int>& density, double minY, double maxY,
                        double dy) {
    int yind;
    for (auto &atom:_atoms) {
        yind = floor((atom.get_position()[1]-minY)/dy);
        if (yind>0 && yind < velocity.size()){
            velocity[yind]+=atom.get_velocity()[0];
            density[yind]+=1;
        }
    }


}
bool Frame::readFacet(const char* fileName) {
    std::fstream fp;
    // FILE* fp;
    fp.open(fileName,std::ios::in);
    std::string line;
    // fp = fopen(fileName, "r");
    if (!fp.is_open()) {
        //perror("Could not open dihedral file");
        //:std::cout<<"dihedral file: "<<fileName<<std::endl;
        // printf(stderr, fileName);
        // printf(stderr, "\n");
        return false;
    }
    // skip the first 3 lines
    for (int i = 0; i < 3; i++) getline(fp,line);//fgets(line, 50, fp);  // read line n
    // read the number of bonds
    int numDihedrals;
    getline(fp,line);
    sscanf(line.c_str(), "%d\n", &numDihedrals);
    typedef std::unordered_set<int> inner_type;
    typedef std::unordered_set<inner_type, hash_on_set_of_ints> set_of_unique_sets;
    
    set_of_unique_sets facetSet;
    // std::unordered_set<std::unordered_set<int>,hash_on_sum<int>> facetSet;

    // skip 5
    for (int i = 0; i < 5; i++) getline(fp,line);  // read line n

    int id, type, a, b, c, d;
    for (int i = 0; i < numDihedrals; i++) {
        getline(fp,line);
        sscanf(line.c_str(),"%d %d %d %d %d %d \n", &id, &type, &a, &b, &c, &d);
        //index by 0
        a--;
        b--;
        c--;
        d--;

        if (a<0 || b<0 || c<0 ||d<0){
            std::cout<<"bad dihedral file found dihedral "<< a <<", "<<b<<", "<<c<<", "<<d<<" Next line is:"<<std::endl;
            std::cout<<line;

            throw "bad id number read dihedral";
        }
        
        std::unordered_set<int>
            face1 = {a, b, c},
            face2 = {b, c, d};
        facetSet.insert(face1);
        facetSet.insert(face2);
    }
    // now that we have a set of faces convert it to a vector for fast look up
    // _facets.resize(facetSet.size());
    for (auto &face: facetSet) {
        std::array<int,3> facetArray;
        std::copy(face.begin(),face.end(),facetArray.begin());
        _facets.push_back(facetArray);
    }
    // std::copy(facetSet.begin(), facetSet.end(), _facets.begin());
    fp.close();
    return true;
}
double Frame::area() {
    return Geometry::area(_atoms,_facets);
}

double Frame::radius() {
    int type = 2;
    return Geometry::radius(_atoms,type);
}

