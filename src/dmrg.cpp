#include "utility.h"
#include <tuple>
#include "dmrg.h"
#include <cxxopts.hpp>
#include <chrono>
int main(int argc, char** argv){
    cxxopts::Options options("MyProgram", "One line description of MyProgram");
    options.add_options()
        ("s,site", "Site types", cxxopts::value<std::string>())
        ("b,bc", "Boundary condition types", cxxopts::value<std::string>())
        ("n", "# of sites", cxxopts::value<int>())
        ("d", "Max dimension", cxxopts::value<int>())
        ("c,cutoff", "Cutoff", cxxopts::value<float>())
        ("t,tol", "Tolerance", cxxopts::value<float>())
        ("theta", "theta", cxxopts::value<float>())
        ("u,penalty", "Penalty size", cxxopts::value<float>())
        ("nstates", "# of states to solve for", cxxopts::value<int>())
        ("analysis", "analysis?", cxxopts::value<int>())
            ;
    auto result = options.parse(argc, argv);
    //fixme: print parameters used
    auto params = std::make_tuple(
            result["s"].as<std::string>(),
            result["b"].as<std::string>(),
            result["n"].as<int>(),
            result["d"].as<int>(),
            result["c"].as<float>(),
            result["t"].as<float>(),
            result["theta"].as<float>(),
            result["u"].as<float>(),
            result["nstates"].as<int>(),
            result["analysis"].as<int>());

    // fixme: disentangle case from task

    if(std::get<0>(params) == "golden"){
        auto dmrg_ = DMRG<Golden>(params);
        if (std::get<9>(params) == 1){
            dmrg_.analyze();
        } else {
            auto start = std::chrono::high_resolution_clock::now();
            dmrg_.run();
            auto stop = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;
        }
        return 0;
    }else if(std::get<0>(params) == "haagerup"){
        auto dmrg_ = DMRG<Haagerup>(params);
        if (std::get<9>(params) == 1) {
            dmrg_.analyze();
        } else {
            dmrg_.run();
        }
        return 0;
    }else if(std::get<0>(params) == "haagerupQN"){
        auto dmrg_ = DMRG<HaagerupQ>(params);
        if (std::get<9>(params) == 1) {
            dmrg_.analyze();
        } else {
            dmrg_.run();
        }
        return 0;
    }else {
        throw std::invalid_argument("Invalid site type.");
    }
}