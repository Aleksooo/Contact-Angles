#include <iostream>
#include <filesystem>
#include "handler.hpp"
#include "struct/Atom.hpp"
#include "struct/Molecule.hpp"
#include "struct/System.hpp"
#include "IO/io_gro.hpp"
#include "shapes/Box.hpp"
#include "shapes/Cylinder.hpp"
#include "shapes/AntiCylinder.hpp"
#include "alg/generator.hpp"
#include "random"

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());

    System sys("Test", vec({3, 3, 3}));

    std::filesystem::path decane_inp("../ff/gromos/gro/decane.gro");
    Molecule decane = read_mol(decane_inp);
    decane.set_atoms_to_center();

    Cylinder cylinder(sys.get_center(), 1, 3, vec({1, 0, 0}));

    std::cout << "Inserting decane..." << std::endl;
    for (size_t mol_id = 1; mol_id <= 26; mol_id++) {
        // std::cout << mol_id << std::endl;
        insert_mol_into_shape(sys, cylinder, decane, mol_id, gen);
    }


    std::filesystem::path water_inp("../ff/gromos/gro/water.gro");
    Molecule water = read_mol(water_inp);
    water.set_atoms_to_center();

    // AntiCylinder anticylinder(vec({1.5, 1.5, 1.5}), 1, 3, vec({1, 0, 0}), vec({1.5, 1.5, 1.5}), sys.box);
    AntiCylinder anticylinder(cylinder, Box(sys.get_center(), sys.box));

    // std::cout << anticylinder.borders_center << std::endl;
    // std::cout << anticylinder.borders << std::endl;

    std::cout << "Inserting water..." << std::endl;
    for (size_t mol_id = 1; mol_id <= 580; mol_id++) {
        // std::cout << mol_id << std::endl;
        insert_mol_into_shape(sys, anticylinder, water, mol_id, gen);
    }

    push_atoms_apart(sys, gen);

    std::filesystem::path outp("../decane_water_cylinder_test.gro");
    write_sys(sys, outp);

    return 0;
}
