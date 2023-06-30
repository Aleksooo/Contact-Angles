#include <iostream>
#include <vector>
#include <numeric>
#include "Atom.hpp"

#include "../handler.hpp"

#ifndef MOLECULE
#define MOLECULE

class Molecule {
private:
    int id = 1;
    std::string name = "MOL";
    std::vector<Atom> atoms;

public:
    // RAII
    Molecule();
    Molecule(const Molecule& lhs);
    Molecule& operator=(const Molecule& lhs);
    Molecule(Molecule&& rhs);
    Molecule& operator=(Molecule&& rhs);

    // Methods
    vec get_center();
    void add_atom(const Atom& lhs);
    std::vector<Atom> get_atoms();
    void set_atoms_to_center();
    void set_id(int id_);
    void set_name(std::string name_);
};

#endif
