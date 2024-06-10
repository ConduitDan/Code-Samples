
#include "lammpsPostProcessor.h"
#include "observers.h"
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

void write_one_min_movie(std::string projectFolder);
std::tuple<int,int,int> findTS(std::string projectFolder);
int main(int argc , char**argv) {

    if (argc<2){
        std::cout<<"Wrong number of arguments ("<<argc<<") please provide the project folder";
    }
    std::string projectFolder = argv[1];
    LammpsPostProcessor LPP(projectFolder);
    std::vector<FrameObserver*> obvs;
    std::tuple<int,int,int> TSData = findTS(projectFolder);
    int initTS = std::get<0>(TSData);
    int finalTS = std::get<1>(TSData);
    int incrementTS = std::get<2>(TSData);
    std::cout<<"Found First Dump at "<<initTS<<" last dump at "<<finalTS<<" with increments of "<<incrementTS<<std::endl;


    if (argc>=3){
        std::string command = argv[2];

        if (command == "flipTest") {
        
            LPP.read_data(0,0,1);
            if (LPP.check_flips()) {
                std::cout<< "Test Passed"<<std::endl;
            } else std::cout <<"FAILED"<<std::endl;
            return 1;
        } else if (command == "newShell") {
            std::tuple<int,int,int> TSData = findTS(projectFolder);
            int initTS = std::get<0>(TSData);
            int finalTS = std::get<1>(TSData);
            int incrementTS = std::get<2>(TSData);
            std::cout<<"Found First Dump at "<<initTS<<" last dump at "<<finalTS<<" with increments of "<<incrementTS<<std::endl;
            LPP.read_data(finalTS,finalTS,incrementTS);
            LPP.write_frame_to_input();

        } else if (command == "image") {

            TSData=std::make_tuple(std::get<1>(TSData),std::get<1>(TSData),std::get<2>(TSData));
            obvs.push_back(new MovieObserver(projectFolder));
        } else if (command == "visc_movie") {
            obvs.push_back(new ViscMovieObserver(projectFolder,false));
            obvs.push_back(new VelocityProfileObserver(projectFolder,(std::get<0>(TSData)-std::get<1>(TSData))/(std::get<2>(TSData))/2));

        } else if (command == "visc") {
            // obvs.push_back(new ViscMovieObserver(projectFolder,false));
            obvs.push_back(new VelocityProfileObserver(projectFolder,(std::get<0>(TSData)-std::get<1>(TSData))/(std::get<2>(TSData))/2));

        } else if (command == "test") {
            double k = atof(argv[3]);
            double lmin = atof(argv[4]);
            double lmax = atof(argv[5]);
            int r = atoi(argv[6]);
            obvs.push_back(new ForceObserver("power",projectFolder,k,lmin,lmax,r));

        } else if (command == "test_ph") {
            double k = atof(argv[3]);
            double lmin = atof(argv[4]);
            double lmax = atof(argv[5]);
            int r = atoi(argv[6]);
            obvs.push_back(new ForceObserver("power_harmonic",projectFolder,k,lmin,lmax,r));
        } else if (command == "imageVel") {
            TSData=std::make_tuple(std::get<1>(TSData),std::get<1>(TSData),std::get<2>(TSData));
            obvs.push_back(new ViscMovieObserver(projectFolder,false));
        } else if (command == "one_min_movie") {
            write_one_min_movie(projectFolder);
            return 0;
        }
    } else {
        obvs.push_back(new MovieObserver(projectFolder));
        obvs.push_back(new AreaObserver(projectFolder));
        // obvs.push_back(new VolumeObserver(projectFolder));
        obvs.push_back(new RadiusObserver(projectFolder));
        obvs.push_back(new SwapsObserver(projectFolder));
    }
    try {
    LPP.post_process_frames(TSData,obvs);
    } catch (const char * str) {
        std::cout << "Exception: " << str << std::endl;
    }
    for (int i = 0; i<obvs.size(); i++){
	std::cout<<"Finishing Up Observers"<<std::endl;
	obvs[i]->finish();
	delete obvs[i];
    }
    return 0;

}

std::tuple<int,int,int> findTS(std::string projectFolder){
    std::vector<int> fileNumbers;
    for (const auto & entry : fs::directory_iterator(projectFolder+"/trajectory/")) {
        int temp;
	fs::path fn = entry.path().filename();
        sscanf(fn.c_str(),"%d.dump",&temp);
        fileNumbers.push_back(temp);
    }
    std::tuple<int,int,int> out;
    sort(fileNumbers.begin(), fileNumbers.end());
    out = {fileNumbers[0],fileNumbers[fileNumbers.size()-1],fileNumbers[1]-fileNumbers[0]};
    return out;
    
}
void write_one_min_movie(std::string projectFolder){
    std::vector<int> fileNumbers;
    for (const auto & entry : fs::directory_iterator(projectFolder + "/images/")) {
        int temp;
        fs::path fn = entry.path().filename();
        sscanf(fn.c_str(),"%d.png",&temp);
        fileNumbers.push_back(temp);
	
    }
    sort(fileNumbers.begin(), fileNumbers.end());
    int end = fileNumbers[fileNumbers.size()-1]; 
    int start = 0;
    if (end-1800>0) start = end -1800;
    std::fstream imageFile;
    imageFile.open(projectFolder+"/frame_file.txt",std::fstream::out);

    for (int i = start; i<=end; i++){
         imageFile<<i<<".png"<<std::endl;
    }
    imageFile.close();
   // system("ffmpeg -i frame_file.txt movie.mp4");

}

