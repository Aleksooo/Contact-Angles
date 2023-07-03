#include "System.hpp"

System::System() {}

System::System(const System& lhs) : title(lhs.title), box(lhs.box), atoms(lhs.atoms) {}

System& System::operator=(const System& lhs) {
    System t(lhs);
    std::swap(title, t.title);
    std::swap(box, t.box);
    std::swap(atoms, t.atoms);

    return *this;
}

System::System(System&& rhs) : title(rhs.title), box(rhs.box), atoms(rhs.atoms) {
    rhs.title = "System";
    rhs.box = vec();
    rhs.atoms.clear();
}

System& System::operator=(System&& rhs) {
    System t(std::move(rhs));
    std::swap(title, t.title);
    std::swap(box, t.box);
    std::swap(atoms, t.atoms);

    return *this;
}


void System::add_mol(const Molecule& mol_) {
    for (Atom a : mol_.atoms) {
        atoms.push_back(a);
    }
}
