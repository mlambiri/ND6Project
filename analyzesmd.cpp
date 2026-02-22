// C++11
// analyze steered molecular dynamics trajectories
// written by Jack 20211103
// edited 20211202
// compile it like this for $$$$$:
// g++ -I /u/project/ana/jtfuller/Programs/eigen-3.3.9/ -std=c++11 -O3 analyzesmd.cpp -o analyzesmd

#include <eigen3/Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <sys/stat.h>  // for stat (counting directories)
#include <unistd.h>  // for getopt

const double gas_constant = 8.3144626 / 4184;  // kcal/K/mol


int main( int argc, char * argv[]) {

    // argument parsing
    std::string jobname;
    int maxtraj = 0;
    double rstepsize = 0.1;
    double temperature = 310;
    int a;
    opterr = 0;
    while ( (a = getopt(argc, argv, "m:s:T:")) != -1 ) {
        switch (a) {
            case 'm':
                maxtraj = std::stoi(optarg);
                break;
            case 's':
                rstepsize = std::stod(optarg);
                break;
            case 'T':
                temperature = std::stod(optarg);
                break;
            default:
                throw a;  // int
        }
    }
    if (optind < argc) {
        jobname = argv[optind];
    } else {
        throw argc;  // int
    }

    // counting traj directories
    if (maxtraj == 0) {
        struct stat info;
        while( !stat(("traj" + std::to_string(maxtraj + 1)).c_str(), &info) ) {
            maxtraj++;
        }
    }

    // counting completed trajectories
    int numtraj = 0;
    for (int trajid=1; trajid <= maxtraj; trajid++) {
        struct stat info;
        if (!stat(("traj" + std::to_string(trajid) + "/" + jobname + "_traj" + std::to_string(trajid)
                + ".colvars.traj").c_str(), &info)) {
            numtraj++;
        }
    }

    // read colvars conf file
    int cvfreq;
    double rbegin;
    double rend;
    int numsteps;
    double force_constant;
    std::ifstream colvarsconffile(("traj1/" + jobname + "_traj1.colvars.conf").c_str());
    if (colvarsconffile.is_open()) {
        std::string line;
        while (std::getline(colvarsconffile, line)) {
            if (line.size() > 0) {
                line = line.substr(line.find_first_not_of(' '));  // remove preceding whitespace
                line = line.substr(0, line.find_first_of('#'));  // remove comments
                if (line.substr(0, 20) == "colvarstrajfrequency") {
                    cvfreq = std::stoi(line.substr(21));
                } else if (line.substr(0, 13) == "forceconstant") {
                    force_constant = std::stod(line.substr(14));
                } else if (line.substr(0, 7) == "centers") {
                    rbegin = std::stod(line.substr(8));
                } else if (line.substr(0, 13) == "targetcenters") {
                    rend = std::stod(line.substr(14));
                } else if (line.substr(0, 14) == "targetnumsteps") {
                    numsteps = std::stoi(line.substr(15));
                }
            }
        }
        colvarsconffile.close();
    } else {
        throw colvarsconffile;  // std::ifstream
    }

    // read colvars traj files
    int timesteps = numsteps / cvfreq + 1;
    Eigen::ArrayXXd rlists (numtraj, timesteps);
    Eigen::ArrayXd centerlist (timesteps);
    Eigen::ArrayXXd expworklists (numtraj, timesteps);
    {
        int trajnum = 0;
        for (int trajid=1; trajid <= maxtraj; trajid++) {
            std::ifstream trajfile(("traj" + std::to_string(trajid) + "/" + jobname + "_traj"
                                    + std::to_string(trajid) + ".colvars.traj").c_str());
            if (trajfile.is_open()) {
                std::string line;
                int timestep = 0;
                while (std::getline(trajfile, line)) {
                    if (line.substr(0, 1) != "#") {
                        rlists(trajnum, timestep) = std::stod(line.substr(12, 24));
                        if (trajnum == 0) {
                            centerlist(timestep) = std::stod(line.substr(36, 23));
                        }
                        expworklists(trajnum, timestep) = std::stod(line.substr(59, 22));
                        timestep++;
                    }
                }
                if (timestep != timesteps) {
                    std::cout << " Problem encountered with incomplete trajectory in file traj" + std::to_string(trajid)
                            + "/" + jobname + "_traj" + std::to_string(trajid) + ".colvars.traj" << std::endl;
                    throw timestep;
                }
                trajfile.close();
                trajnum++;
            } else {
                std::cout << " Problem encountered opening file traj" + std::to_string(trajid) + "/" + jobname +
                             "_traj" + std::to_string(trajid) + ".colvars.traj\n Skipping..." << std::endl;
            }
        }
        if (trajnum != numtraj) {
            std::cout << " Problem with number of trajectories" << std::endl;
            throw trajnum;
        }
    }

    // calculating potential of mean force
    int rsteps = (rend - rbegin) / rstepsize + 1;
    expworklists = exp( -expworklists / gas_constant / temperature );
    Eigen::ArrayXd expworksum = expworklists.colwise().sum();
    Eigen::ArrayXd numerator = Eigen::ArrayXd::Zero(rsteps);
    Eigen::ArrayXd denominator = Eigen::ArrayXd::Zero(rsteps);
    for (int rnum=0; rnum < rsteps; rnum++) {
        double r = rbegin + rnum * rstepsize;
        for ( int timestep=0; timestep < timesteps; timestep++) {
            double expv = exp( -force_constant / 2 * pow(r - centerlist(timestep), 2.0) / gas_constant / temperature );
            denominator(rnum) += expv / expworksum(timestep);
            for (int trajnum=0; trajnum < numtraj; trajnum++) {
                if ( std::abs(rlists(trajnum, timestep) - r) < rstepsize / 2 ) {
                    numerator(rnum) += expworklists(trajnum, timestep) / expworksum(timestep);
                }
            }
        }
    }
    Eigen::ArrayXd pmf = -gas_constant * temperature * log( numerator / denominator / numtraj );

    // write pmf to file
    std::ofstream pmffile ( (jobname + ".pmf").c_str() );
    if (pmffile.is_open()) {
        pmffile << "# " << jobname << ' ' << numtraj << "-trajectory pmf" << std::endl;
        pmffile << "# r (Å), G (kcal/mol)" << std::endl;
        for (int rnum=0; rnum < rsteps; rnum++) {
            pmffile << rbegin + rnum * rstepsize << ',' << pmf(rnum) << std::endl;
        }
        pmffile.close();
    } else {
        throw pmffile;  // std::ofstream
    }

}
