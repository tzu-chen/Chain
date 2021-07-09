#ifndef CHAIN_PARAMETERS
#define CHAIN_PARAMETERS

#include <filesystem>

// Path to the directory under which subdirectories such as kPath/pgs will be created to store data files.
// Default: /Chain.
std::filesystem::path kPath = std::filesystem::current_path().parent_path();

#endif