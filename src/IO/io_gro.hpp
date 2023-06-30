#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <algorithm>
#include <array>
#include <vector>
#include "../comps/Atom.hpp"
#include "../comps/Molecule.hpp"
#include "../comps/System.hpp"

#include "../handler.hpp"

#ifndef IO
#define IO

Molecule read_mol(std::filesystem::path file);

#endif
