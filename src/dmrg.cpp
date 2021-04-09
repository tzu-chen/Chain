#include <filesystem>
#include "utility.cpp"
#include "dmrg_progress.cpp"
#include "dmrg.h"

int main(int argc, char** argv){
    // fixme: distangle case from task
    if(std::string(argv[1]) == "golden"){
        if(argv[10] != nullptr){
            DMRG<Golden>(argv).analyze();
        } else {
            DMRG<Golden>(argv).run();
        }
        return 0;
    }else if(std::string(argv[1]) == "haagerup"){
        if(argv[10] != nullptr){
            DMRG<HaagerupQ>(argv).analyze();
        } else {
            DMRG<HaagerupQ>(argv).run();
        }
        return 0;
    }
}