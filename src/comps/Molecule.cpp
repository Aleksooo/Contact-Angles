#include "Molecule.hpp"

Molecule::Molecule() {}

Molecule::Molecule(const Molecule& lhs) : id(lhs.id), name(lhs.name), atoms(lhs.atoms) {}

Molecule& Molecule::operator=(const Molecule& lhs) {
    Molecule t(lhs);
    std::swap(id, t.id);
    std::swap(name, t.name);
    std::swap(atoms, t.atoms);

    return *this;
}

Molecule::Molecule(Molecule&& rhs) : id(rhs.id), name(rhs.name), atoms(rhs.atoms) {
    rhs.id = 1;
    rhs.name = "MOL";
    rhs.atoms.clear();
}

Molecule& Molecule::operator=(Molecule&& rhs) {
    Molecule t(std::move(rhs));
    std::swap(id, t.id);
    std::swap(name, t.name);
    std::swap(atoms, t.atoms);

    return *this;
}


vec Molecule::get_center() {
    vec sum = std::accumulate(atoms.begin(), atoms.end(), vec(), [](vec a, Atom b){ return a + b.xyz; });
    return sum / atoms.size();
}

void Molecule::set_atoms_to_center() {
    vec center = get_center();
    for (Atom& atom : atoms) {
        atom.xyz = (atom.xyz - center);
    }
}
