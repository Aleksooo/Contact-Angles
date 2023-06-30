#include "io_gro.hpp"

Molecule read_mol(std::filesystem::path rfile) {
    std::ifstream file(rfile);
    Molecule mol;
    Atom atom;

    std::string line, subline;
    std::vector<std::string> args;
    std::array<int, 7> mask = {5, 5, 5, 5, 8, 8, 8};
    size_t counter = 0;
    std::getline(file, line); // skip title
    std::getline(file, line);
    size_t size = std::stoul(line);
    std::cout << size << std::endl;

    for (size_t i = 0; i < size; i++) {
        std::getline(file, line);
        std::cout << line << std::endl;
        for (size_t p = 0; p < line.size(); p += mask[counter]) {
            // std::cout << line.substr(p, mask[counter]) << std::endl;
            args.push_back(std::string(line.substr(p, mask[counter])));
            args[counter].erase(std::remove_if(args[counter].begin(), args[counter].end(), isspace), args[counter].end());
            std::cout << args[counter] << std::endl;
            counter++;
        }

        mol.set_id(std::stoi(args[0]));
        mol.set_name(args[1]);
        atom.set_name(args[2]);
        atom.set_id(std::stoi(args[3]));
        atom.set_XYZ(vec({std::stod(args[4]), std::stod(args[5]), std::stod(args[6])}));

        mol.add_atom(atom);
        args.clear();
        counter = 0;
    }

    return mol;

}
