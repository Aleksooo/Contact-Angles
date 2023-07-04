#include "generator.hpp"

bool check_simple_intersections(System& sys, vec point, double min_dist2) {
    for (const Atom& a : sys.atoms) {
        if (dist2_pbc(a.xyz, point, sys.box) < min_dist2) {
            return true;
        }
    }

    return false;
}

void insert_mol_into_shape(
    System& sys,
    const Shape& shape,
    Molecule mol,
    size_t mol_id,
    std::mt19937& gen,
    int insertion_limit,
    int rotation_limit,
    double min_dist2,
    double package
) {
    bool overlap = true;
    int insertion_counter = 0;
    int rotation_counter = 0;

    vec new_point;

    // Try to insert
    while (overlap && (insertion_counter < insertion_limit)) {
        insertion_counter++;

        new_point = shape.generate_point(gen);
        overlap = check_simple_intersections(
            sys,
            new_point,
            (2 - insertion_counter / insertion_limit) * package * mol.size
        );

        // rotating
        while (overlap && (rotation_counter < rotation_limit)) {
            rotation_counter++;
            random_mol_rotation(mol, gen);

            overlap = check_simple_intersections(sys, new_point, min_dist2);
        }
        rotation_counter = 0;
    }

    if (insertion_counter == insertion_limit) { std::cerr << "FAIL!" << std::endl; }

    mol.set_atoms_to_point(new_point);
    sys.add_mol(mol, mol_id);
}


void push_atoms_apart(
    System& sys,
    std::mt19937& gen,
    double min_dist2,
    double max_dist2,
    int iteration_lim
) {
    double min_dist = sqrt(min_dist2);
    double max_dist = sqrt(max_dist2);
    std::uniform_real_distribution<double> dist(min_dist, max_dist);

    int iter = 0;
    bool overlap = true;
    int overlap_counter = 0;

    while (overlap && (iter < iteration_lim)) {
        iter++;
        overlap = false;
        overlap_counter = 0;

        for (size_t i = 0; i < sys.atoms.size() - 1; i++) {
            for (size_t j = i + 1; j < sys.atoms.size(); j++) {
                if (sys.atoms[i].mol_id != sys.atoms[j].mol_id) {
                    vec delta = delta_pbc(sys.atoms[i].xyz, sys.atoms[j].xyz, sys.box);

                    if (norm(delta) < min_dist) {
                        sys.atoms[i].xyz = sys.atoms[i].xyz - delta * dist(gen) / norm(delta);
                        overlap = true;
                        overlap_counter++;
                    }
                }
            }
        }
    }

    std::cout << overlap_counter << " overlaps detected" << std::endl;
    sys.apply_pbc();
}


void random_mol_rotation(Molecule& mol, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(0, M_PI);
    double psi = 2 * dist(gen);
    double cos_psi = std::cos(psi);
    double sin_psi = std::sin(psi);

    double theta = dist(gen);
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    double phi = 2 * dist(gen);
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    // Formula from http://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf
    for (Atom& a : mol.atoms) {
        a.xyz = vec({
            cos_theta * cos_phi * a.xyz[0] + (sin_psi * sin_theta * cos_phi - cos_psi * sin_phi) * a.xyz[1] + (cos_psi * sin_theta * cos_phi + sin_psi * sin_phi) * a.xyz[2],
            cos_theta * sin_phi * a.xyz[0] + (sin_psi * sin_theta * sin_phi + cos_psi * cos_phi) * a.xyz[1] + (cos_psi * sin_theta * sin_phi - sin_psi * cos_phi) * a.xyz[2],
            -sin_theta * a.xyz[0] + sin_psi * cos_theta * a.xyz[1] + cos_psi * cos_theta * a.xyz[2]
        });
    }
}
