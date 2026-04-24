#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "hypergraph_motif.hpp"


int main(int argc, char* argv[])
{
    if (argc < 2) { std::cerr << "Error: missing input file.\n"; return 1;}

    std::string input {argv[1]};

    std::string output {};
    int motif_size {4};

    for (int i = 2; i < argc; i += 2){
        std::string arg {argv[i]};

        if (arg == "-o") {
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
        output = nameWithoutExt + "_result.txt";
    }

    Hypergraph H {read_hyperedge_list(input)};

    std::cout << "Nodes: " << H.node_count() << ". Edges: " << H.edge_count() << ".\n";

    auto size_hist {H.edge_size_count()};
    for (std::size_t i {0}; i < std::size(size_hist); ++i){
        std::cout << "Size " << i+2 << ":\t" << size_hist[i] << '\n';
    }

    std::cout << "\nComputing frequencies of order " << motif_size << " motifs...\n\n";

    std::ofstream outfile {output};
    if (!outfile) { std::cerr << "Error opening output file.\n"; return 1; }

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::_V2::system_clock::time_point end {};

    if (motif_size == 3) {
        auto counts {H.motif_count_3()};
        end = std::chrono::high_resolution_clock::now();
        for (std::size_t i {}; i < std::size(counts); ++i){
            outfile << "Type " << i+1 << ":\t" << counts[i] << '\n';
        }
    }

    else {
        auto counts {H.motif_count_4()};
        end = std::chrono::high_resolution_clock::now();
        for (std::size_t i {}; i < std::size(counts); ++i){
            outfile << "Type " << i+1 << ":\t" << counts[i] << '\n';
        }
    }

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Time taken: " << duration << " ms\n";
    std::cout << "Results saved under " << output << '\n';

    return 0;
}