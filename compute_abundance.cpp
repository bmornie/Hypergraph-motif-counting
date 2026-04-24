#include <iostream>
#include <fstream>
#include <string>

#include "hypergraph_motif.hpp"
#include "random_models.hpp"


int main(int argc, char* argv[])
{
    if (argc < 2) { std::cerr << "Error: missing input file.\n"; return 1;}

    std::string input {argv[1]};

    std::string model {"conf"};
    int samples {100};
    double epsilon {4.0};
    int swaps {-1};
    std::string output {};
    int motif_size {4};

    for (int i = 2; i < argc; i += 2){
        std::string arg {argv[i]};

        if (arg == "-m") {
            model = argv[i+1];
            if (model != "sp" && model != "conf") { std::cerr << "Error: model should be either 'sp' (size-preserving) or 'conf' (configuration).\n"; return 1; }
        }
        else if (arg == "-n") {
            samples = {std::stoi(argv[i+1])};
            if (samples <= 0) { std::cerr << "Error: number of samples should be a positive integer.\n"; return 1; }
        }
        else if (arg == "-e") {
            epsilon = {std::stod(argv[i+1])};
        }
        else if (arg == "-t") {
            swaps = {std::stoi(argv[i+1])};
            if (swaps <= 0) { std::cerr << "Error: number of hyperedge swaps should be a positive integer.\n"; return 1; }
        }
        else if (arg == "-o") {
            output = argv[i+1];
        }
        else if (arg == "-s") {
            int s {std::stoi(argv[i+1])};
            if (s != 3 && s != 4) { std::cerr << "Error: motif size can only be 3 or 4.\n"; return 1; }
            motif_size = s;
        }
    }

    if (output.empty()) {
        std::string::size_type slashPos = input.find_last_of("/\\");
        std::string filename = (slashPos == std::string::npos) ? input : input.substr(slashPos + 1);
        std::string::size_type dotPos = filename.find_last_of('.');
        std::string nameWithoutExt = (dotPos == std::string::npos) ? filename : filename.substr(0, dotPos);
        output = nameWithoutExt + "_abundance_" + model + ".txt";
    }

    Hypergraph H {read_hyperedge_list(input)};

    std::cout << "Nodes: " << H.node_count() << ". Edges: " << H.edge_count() << ".\n";

    auto size_hist {H.edge_size_count()};
    for (std::size_t i {0}; i < std::size(size_hist); ++i){
        std::cout << "Size " << i+2 << ":\t" << size_hist[i] << '\n';
    }

    std::string modelname {model == "conf" ? "configuration" : "size-preserving"};
    std::cout << "\nComputing abundance based on motifs of order " << motif_size << " under " << modelname << " model...\n\n";

    std::ofstream outfile {output};
    if (!outfile) { std::cerr << "Error opening output file.\n"; return 1; }

    int m {model == "sp" ? 0 : 1};
    auto abundance {motif_abundance(H, motif_size, m, samples, swaps, epsilon)};

    if (motif_size == 3) {
        auto counts {H.motif_count_3()};

        if (counts.size() != abundance.size()) { std::cerr << "Something went wrong...\n"; return 1; }

        for (size_t i {}; i < counts.size(); ++i) {
            outfile << "Type " << i+1 << ":\t" << counts[i] << '\t' << abundance[i] << '\n';
        }
    }

    else {
        auto counts {H.motif_count_4()};

        if (counts.size() != abundance.size()) { std::cerr << "Something went wrong...\n"; return 1; }

        for (size_t i {}; i < counts.size(); ++i) {
            outfile << "Type " << i+1 << ":\t" << counts[i] << '\t' << abundance[i] << '\n';
        }
    }

    std::cout << "Results saved under " << output << '\n';

    return 0;
}