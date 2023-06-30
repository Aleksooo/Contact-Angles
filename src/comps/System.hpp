#include <iostream>
#include <vector>
#include "Atom.hpp"

#include "../handler.hpp"

#ifndef SYSTEM
#define SYSTEM

class System {
private:
    std::string title = "System";
    vec box;
    std::vector<Atom> atoms;

public:
    // RAII
    System();
    System(const System& lhs);
    System& operator=(const System& lhs);
    System(System&& rhs);
    System& operator=(System&& rhs);

    // Methods
};

#endif
