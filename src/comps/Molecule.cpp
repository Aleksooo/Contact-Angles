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
    vec sum = std::accumulate(atoms.begin(), atoms.end(), vec(), [](auto a, auto b){ return a + b.get_XYZ(); });
    return sum / atoms.size();
}

void Molecule::add_atom(const Atom& lhs) {
    atoms.push_back(Atom(lhs));
}

std::vector<Atom> Molecule::get_atoms() {
    return atoms;
}

void Molecule::set_atoms_to_center() {
    vec center = get_center();
    for (auto& atom : atoms) {
        atom.set_XYZ(atom.get_XYZ() - center);
    }
}

void Molecule::set_id(int id_) {
    id = id_;
}

void Molecule::set_name(std::string name_) {
    name = name_;
}
