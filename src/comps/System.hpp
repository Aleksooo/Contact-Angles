#include <iostream>
#include <vector>
#include "Atom.hpp"
#include "Molecule.hpp"

#include "../handler.hpp"

#ifndef SYSTEM
#define SYSTEM

class System {
public:
    std::string title = "System";
    vec box;
    std::vector<Atom> atoms;

    // RAII
    System();
    System(const System& lhs);
    System& operator=(const System& lhs);
    System(System&& rhs);
    System& operator=(System&& rhs);

    // Methods
    void add_mol(const Molecule& mol_);
};

#endif
