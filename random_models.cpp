#include "random_models.hpp"

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <cmath>

#include "hypergraph_motif.hpp"
#include "utilities/random.hpp"


Hypergraph sample_sizepreserving(const Hypergraph& H_orig)
{
    auto size_count {H_orig.edge_size_count()};
    ID n {static_cast<ID>(H_orig.node_count())};

    Hypergraph H {};

    for (size_t s {}; s < 3; ++s) {
        size_t size {s+2};
        auto sc {size_count[s]};

        size_t t = 0;
        while (t < sc) {
            std::vector<ID> edge {rand_int(static_cast<ID>(0), n-1, size)};
            if (!H.has_edge(edge)) {
                H.add_edge(edge);
                ++t;
            }
        }
    }
    return H;
}


static void _pairwise_shuffle(std::vector<ID>& edge1, std::vector<ID>& edge2)
{
    std::vector<ID> new1;
    std::vector<ID> new2;
    new1.reserve(edge1.size());
    new2.reserve(edge2.size());

    std::vector<ID> unique_nodes {};
    for (auto n : edge1){
        if (std::find(edge2.begin(), edge2.end(), n) == edge2.end()){
            unique_nodes.push_back(n);
        }
        else {
            new1.push_back(n);
            new2.push_back(n);
        }
    }
    for (auto n : edge2){
        if (std::find(edge1.begin(), edge1.end(), n) == edge1.end())
            unique_nodes.push_back(n);
    }

    while (unique_nodes.size() > 0){
        ID n {pop_random(unique_nodes)};
        if (new1.size() < edge1.size())
            new1.push_back(n);
        else
            new2.push_back(n);
    }

    std::sort(new1.begin(), new1.end());
    std::sort(new2.begin(), new2.end());

    edge1 = std::move(new1);
    edge2 = std::move(new2);
}

Hypergraph sample_configuration(const Hypergraph& H_orig, int swaps)
{
    std::vector<std::vector<ID>> edge_list {H_orig.edge_list()};
    std::unordered_set<std::vector<ID>, VectorHash> edge_set {edge_list.begin(), edge_list.end()};
    auto m {edge_list.size()-1};

    int t {};
    while (t < swaps){
        auto i {rand_int<std::size_t>(0, m)};
        auto j {rand_int<std::size_t>(0, m)};
        while (i == j)
            j = rand_int<std::size_t>(0, m);

        auto new_edge1 {edge_list[i]};
        auto new_edge2 {edge_list[j]};
        
        _pairwise_shuffle(new_edge1, new_edge2);

        if (!edge_set.contains(new_edge1) && !edge_set.contains(new_edge2)){
            edge_set.erase(edge_list[i]);
            edge_set.erase(edge_list[j]);
            edge_set.insert(new_edge1);
            edge_set.insert(new_edge2);
            edge_list[i] = std::move(new_edge1);
            edge_list[j] = std::move(new_edge2);
            ++t;
        }
    }

    Hypergraph H {};
    for (const auto& edge : edge_list)
        H.add_edge(edge);
    
    return H;
}


static std::vector<double> _motif_3_abundance(const Hypergraph& H, int model, int samples, int swaps, double epsilon)
{ 
    auto motif_counts {H.motif_count_3()};
    
    std::vector<long long> sum (6);

    for (int s {}; s < samples; ++s){
        Hypergraph Hs {model ? sample_configuration(H, swaps) : sample_sizepreserving(H)};
        auto counts {Hs.motif_count_3()};

        for (std::size_t i {}; i < 6; ++i) { sum[i] += counts[i]; }
    }

    std::vector<double> abundance (6);

    for (std::size_t i {}; i < 6; ++i){
        double mean {static_cast<double>(sum[i])/static_cast<double>(samples)};
        abundance[i] = (static_cast<double>(motif_counts[i]) - mean) / (static_cast<double>(motif_counts[i]) + mean + epsilon);
    }

    return abundance;
}

static std::vector<double> _motif_4_abundance(const Hypergraph& H, int model, int samples, int swaps, double epsilon)
{ 
    auto motif_counts {H.motif_count_4()};
    
    std::vector<long long> sum (171);

    for (int s {}; s < samples; ++s){
        Hypergraph Hs {model ? sample_configuration(H, swaps) : sample_sizepreserving(H)};
        auto counts {Hs.motif_count_4()};

        for (std::size_t i {}; i < 171; ++i) { sum[i] += counts[i]; }
    }

    std::vector<double> abundance (171);

    for (std::size_t i {}; i < 171; ++i){
        double mean {static_cast<double>(sum[i])/static_cast<double>(samples)};
        abundance[i] = (static_cast<double>(motif_counts[i]) - mean) / (static_cast<double>(motif_counts[i]) + mean + epsilon);
    }

    return abundance; 
}

std::vector<double> motif_abundance(const Hypergraph& H, int motif_size, int model, int samples, int swaps, double epsilon)
{
    if (model && swaps == -1) { swaps = 10*static_cast<int>(H.edge_count()); }

    switch (motif_size) {
    case 3: return _motif_3_abundance(H, model, samples, swaps, epsilon);
    case 4: return _motif_4_abundance(H, model, samples, swaps, epsilon);
    default: return {0};
    }
}