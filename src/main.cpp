#include <iostream>
#include <filesystem>
#include "handler.hpp"
#include "comps/Atom.hpp"
#include "comps/Molecule.hpp"
#include "IO/io_gro.hpp"
#include "random"

int main() {
    // std::random_device rd;
    // std::mt19937 g(rd());
    // std::uniform_real_distribution<double> dist(-2.0, 2.0);
    // auto gen = [&dist, &g](){ return dist(g); };

    // vec xyz;

    // Molecule mol;
    // Atom atom;
    // for (size_t i = 0; i < 5; i++) {
    //     std::generate(xyz.begin(), xyz.end(), gen);
    //     atom.set_XYZ(xyz);

    //     mol.add_atom(atom);
    // }

    // for (auto a : mol.get_atoms()) {
    //     std::cout << a.get_XYZ() << std::endl;
    // }

    std::filesystem::path path("../ff/gromos/gro/decane.gro");
    Molecule mol = read_mol(path);

    for (auto a : mol.get_atoms()) {
        std::cout << a.get_XYZ() << std::endl;
    }

    return 0;
}
