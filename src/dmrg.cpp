#include "utility.h"
#include <tuple>
#include "dmrg.h"
#include <cxxopts.hpp>
// Uncomment this and the std::chrono lines to enable timing
#include <chrono>

int main(int argc, char** argv) {
    cxxopts::Options options("Chain", "Simulates DMRG on anyon chains.");
    options.add_options()
        ("s,site", "Site types(golden/haagerup/haagerupq)", cxxopts::value<std::string>())
        ("b,bc", "Boundary condition types(p/o/s/sp)", cxxopts::value<std::string>())
        ("l,length", "# of sites", cxxopts::value<int>())
        ("d", "Max bond dimension", cxxopts::value<int>())
        ("c,cutoff", "SVD cutoff", cxxopts::value<float>())
        ("p,precision", "Energy precision", cxxopts::value<float>())
        ("j,couplings", "Couplings", cxxopts::value<std::string>())
        ("u,penalty", "Penalty size", cxxopts::value<float>())
        ("q,charge_", "Charge", cxxopts::value<int>())
        ("n,nstates", "# of states to solve for", cxxopts::value<int>())
        ("m,mode", "DMRG simulation(0)/measurements and analysis(1,2,3)", cxxopts::value<int>())
        ("h,help", "Print usage")
        ;
    auto result = options.parse(argc, argv);
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    auto params = std::make_tuple(
            result["s"].as<std::string>(),
            result["b"].as<std::string>(),
            result["l"].as<int>(),
            result["d"].as<int>(),
            result["c"].as<float>(),
            result["p"].as<float>(),
            result["j"].as<std::string>(),
            result["u"].as<float>(),
            result["q"].as<int>(),
            result["n"].as<int>(),
            result["m"].as<int>());
    // Currently cannot disentangle case from task because the type of dmrg is specific to each case
    if (std::get<0>(params) == "golden") {
        auto dmrg = DMRG<Golden>(params);
        if (std::get<10>(params) == 1) {
            dmrg.Analyze();
        } else if (std::get<10>(params) == 2) {
            dmrg.AnalyzeWithoutRho();
        } else if (std::get<10>(params) == 3) {
            dmrg.Energies();
        } else if (std::get<10>(params) == 4) {
            dmrg.NormalizeEnergies();
        } else if (std::get<10>(params) == 5) {
            dmrg.Simulate(true);
        } else if (std::get<10>(params) == 6) {
            dmrg.Repair();
        } else {
            dmrg.Simulate();
        }
        return 0;
    } else if (std::get<0>(params) == "haagerup") {
        auto dmrg = DMRG<Haagerup>(params);
        if (std::get<10>(params) == 1) {
            auto start = std::chrono::high_resolution_clock::now();
            dmrg.Analyze();
            auto stop = std::chrono::high_resolution_clock::now();
            printf("\nMeasurements total time spent: %gs\n", std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000);
        } else if (std::get<10>(params) == 2) {
            dmrg.AnalyzeWithoutRho();
        } else if (std::get<10>(params) == 3) {
            dmrg.Energies();
        } else if (std::get<10>(params) == 4) {
            dmrg.NormalizeEnergies();
        } else if (std::get<10>(params) == 5) {
            dmrg.Simulate(true);
        } else if (std::get<10>(params) == 6) {
            dmrg.Repair();
        } else {
            dmrg.Simulate();
        }
        return 0;
    } else if (std::get<0>(params) == "haagerupq") {
        auto dmrg = DMRG<HaagerupQ>(params);
        if (std::get<10>(params) == 1) {
            auto start = std::chrono::high_resolution_clock::now();
            dmrg.Analyze();
            auto stop = std::chrono::high_resolution_clock::now();
            printf("\nMeasurements total time spent: %gs\n", std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000);
        } else if (std::get<10>(params) == 2) {
            dmrg.AnalyzeWithoutRho();
        } else if (std::get<10>(params) == 3) {
            dmrg.Energies();
        } else if (std::get<10>(params) == 4) {
            dmrg.NormalizeEnergies();
        } else if (std::get<10>(params) == 5) {
            dmrg.Simulate(true);
        } else if (std::get<10>(params) == 6) {
            dmrg.Repair();
        } else {
            dmrg.Simulate();
        }
        return 0;
    } else {
        throw std::invalid_argument("Invalid site type.");
    }
}