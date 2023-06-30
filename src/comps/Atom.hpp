#include <iostream>
#include <vector>

#include "../handler.hpp"

#ifndef ATOM
#define ATOM

class Atom {
private:
    int mol_id = 1;
    std::string mol_name = "MOL";
    int id = 1;
    std::string name = "ATOM";
    vec xyz;

public:
    // RAII
    Atom();
    Atom(const Atom& lhs);
    Atom& operator=(const Atom& lhs);
    Atom(Atom&& rhs);
    Atom& operator=(Atom&& rhs);

    // Methods
    vec get_XYZ();
    void set_XYZ(const vec& xyz_);
    void set_id(int id_);
    void set_name(std::string name_);
};

#endif
