#include <random>

#include "../struct/Atom.hpp"
#include "../struct/Molecule.hpp"
#include "../struct/System.hpp"
#include "../shapes/Shape.hpp"
#include "../shapes/Box.hpp"
#include "pbc.hpp"
#include <cmath>

#include "../handler.hpp"

bool check_simple_intersections(System& sys, vec point, double min_dist2);

void insert_mol_into_shape(
    System& sys,
    const Shape& shape,
    Molecule mol,
    size_t mol_id,
    std::mt19937& gen,
    int insertion_limit = int(1e5),
    int rotation_limit = 10,
    double min_dist2 = 0.0064,
    double package = 0.4
);

void random_mol_rotation(Molecule& mol, std::mt19937& gen);
