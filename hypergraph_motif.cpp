#include "hypergraph_motif.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <string>

#include "utilities/unordered_dense.h"
#include "motif_lookup.hpp"


void Hypergraph::add_to_index(ID edge_id)
{
    const std::vector<ID>& edge {m_edges[edge_id]};
    auto s {edge.size()};

    if (s == 2){
        std::array<ID,2> tmp {};
        std::copy(edge.begin(), edge.end(), tmp.begin());
        m_edge2_ids.insert(tmp);
        return;
    }
    if (s == 3){
        std::array<ID,3> tmp {};
        std::copy(edge.begin(), edge.end(), tmp.begin());
        m_edge3_ids.insert(tmp);
        return;
    }
    if (s == 4){
        std::array<ID,4> tmp {};
        std::copy(edge.begin(), edge.end(), tmp.begin());
        m_edge4_ids.insert(tmp);
        return;
    }
}

void Hypergraph::add_node(ID n)
{
    if (!m_nodes.contains(n))
        m_nodes[n] = {};
}

void Hypergraph::add_edge(std::vector<ID> nodes)
{
    if (nodes.size() > 4) { return; }
    
    std::sort(nodes.begin(), nodes.end());
    m_edges.push_back(nodes);
    add_to_index(m_next_id);

    for (ID n : nodes)
        m_nodes[n].push_back(m_next_id);

    ++m_next_id;
}

int Hypergraph::degree(ID n, size_t edge_size) const
{
    if (!edge_size)
        return static_cast<int>(m_nodes.at(n).size());

    int deg {};
    for (auto e : m_nodes.at(n))
    {
        if (m_edges.at(e).size() == edge_size)
            ++deg;
    }
    return deg;
}

auto Hypergraph::neighbors(ID n, size_t edge_size) const
{
    std::vector<ID> nbs {};
    for (auto e : m_nodes.at(n))
    {
        if (!edge_size || edge_size == size(e))
        {
            for (auto nb : m_edges.at(e))
            {
                if (nb != n && std::find(nbs.begin(), nbs.end(), nb) == nbs.end())
                    nbs.push_back(nb);
            }
        }
    }
    return nbs;
}

std::vector<int> Hypergraph::degree_list() const
{
    std::vector<int> degrees {};
    degrees.reserve(m_nodes.size());

    for (const auto& [n, e] : m_nodes) { degrees.push_back(static_cast<int>(e.size())); }
    
    return degrees;
}

std::vector<int> Hypergraph::neighbor_count() const
{
    std::vector<int> nbs {};
    nbs.reserve(m_nodes.size());

    for (const auto& [n, _] : m_nodes) {
        nbs.push_back(static_cast<int>(neighbors(n).size()));
    }

    return nbs;
}

std::vector<size_t> Hypergraph::edge_size_count() const
{
    return std::vector<size_t>{m_edge2_ids.size(), m_edge3_ids.size(), m_edge4_ids.size()};
}


std::array<long long, 6> Hypergraph::motif_count_3() const
{
    long long m22 {};   // 2 edges of size 2 forming a wedge
    long long m222 {};  // 3 edges of size 2 forming a triangle
    long long m3 {};	// edge of size 3
    long long m32 {};   // edge of size 3 including a single edge of size 2
    long long m322 {};  // edge of size 3 including a wedge
    long long m3222 {}; // edge of size 3 including a triangle

    // Start from non-induced counts

    for (const auto &[n1, value] : m_nodes) {

        auto nb1 {neighbors(n1, 2)};
        auto d {nb1.size()};
        m22 += static_cast<signed>(d * (d-1) / 2);

        for (size_t i = 0; i < d; ++i) {
            ID n2 {nb1[i]};
            if (n2 < n1) { continue; }

            for (size_t j = i+1; j < d; ++j) {
                ID n3 {nb1[j]};
                if (n3 < n1) { continue; }

                if (has_edge(n2, n3)) { ++m222; }
            }
        }
    }

    for (const auto& edge : m_edges) {
        if (edge.size() != 3) { continue; }

        ++m3;

        const int e01 = has_edge(edge[0], edge[1]);
        const int e02 = has_edge(edge[0], edge[2]);
        const int e12 = has_edge(edge[1], edge[2]);

        m32  += e01 + e02 + e12;
        m322 += e01*e02 + e01*e12 + e02*e12;
        m3222 += e01*e02*e12;
    }

    // Convert to induced counts

    std::array<long long, 6> induced_counts {m22 - 3 * m222 - m322 + 3 * m3222,
                                                m222 - m3222,
                                                m3 - m32 + m322 - m3222,
                                                m32 - 2 * m322 + 3 * m3222,
                                                m322 - 3 * m3222,
                                                m3222};

    return induced_counts;
}

std::array<long long, 171> Hypergraph::motif_count_4() const 
{
    std::array<long long, 171> counts {};

    Hypergraph H {};
    for (const auto &edge : m_edges){
        if (edge.size() == 2)
            H.add_edge(edge);
    }
    H._motif_count_4_size_2(counts);

    for (const auto &edge : m_edges){
        if (edge.size() == 3)
            H.add_edge(edge);
    }
    H._motif_count_4_size_3(counts);

    _motif_count_4_size_4(counts);

    return counts;
}

void Hypergraph::_motif_count_4_size_2(std::array<long long, 171>& counts) const
{
    long long star {};
    long long path {};
    long long tailed_triangle {};
    long long cycle {};
    long long diamond {};
    long long clique {};

    long long triangle {};

    ankerl::unordered_dense::map<uint64_t, long long> wedges {};

    for (const auto &[n, edges] : m_nodes){
        int d {degree(n)};

        star += d*(d-1)*(d-2)/6;

        auto nbs {neighbors(n)};
        size_t N = nbs.size();
        for (size_t i = 0; i < N; ++i) {
            ID u = nbs[i];
            for (size_t j = i+1; j < N; ++j) {
                ID v = nbs[j];

                uint64_t key {};
                if (u>v) { key = (uint64_t(v) << 32) | u; }
                else { key = (uint64_t(u) << 32) | v; }
                ++wedges[key];

                if (n > u || n > v || !has_edge(u,v)) { continue; }

                for (size_t k = j+1; k < N; ++k) {
                    ID w = nbs[k];
                    if (n > w) { continue; }

                    if (has_edge(u, w) && has_edge(v, w))
                        ++clique;
                }
            }
        }
    }

    for (const auto& edge : m_edges){
        ID u = edge[0];
        ID v = edge[1];
        int d1 {degree(u)};
        int d2 {degree(v)};
        uint64_t key = (uint64_t(u) << 32) | v;
        auto w12 {wedges.contains(key)? wedges[key]: 0};

        triangle += w12;
        path += (d1-1)*(d2-1);
        tailed_triangle += w12*(d1 + d2 - 4);
        diamond += w12*(w12-1)/2;
    }

    for (const auto &[n, w] : wedges)
        cycle += w*(w-1)/2;
    
    triangle /= 3;
    path -= 3*triangle;
    tailed_triangle /= 2;
    cycle /= 2;

    counts[0] = star - tailed_triangle + 2*diamond - 4*clique;
    counts[1] = path - 2*tailed_triangle - 4*cycle + 6*diamond - 12*clique;
    counts[2] = tailed_triangle - 4*diamond + 12*clique;
    counts[3] = cycle - diamond + 3*clique;
    counts[4] = diamond - 6*clique;
    counts[5] = clique;
}

void Hypergraph::_motif_count_4_size_3(std::array<long long, 171>& counts) const
{
    for (const auto &edge : m_edges){
        if (edge.size() != 3)  { continue; }
        
        ID n1 {edge[0]};
        ID n2 {edge[1]};
        ID n3 {edge[2]};

        std::vector<ID> edge_nbs {};

        for (ID n : edge) {
            for (ID adj_e: m_nodes.at(n)) {

                ID n4 {};
                int new_count {};

                for (ID x : m_edges[adj_e]) {
                    if (x != n1 && x != n2 && x != n3) {
                        n4 = x;
                        ++new_count;
                    }
                }

                if (new_count != 1) { continue; }

                if (std::find(edge_nbs.begin(), edge_nbs.end(), n4) == edge_nbs.end()) 
                    { edge_nbs.push_back(n4); }
            }
        }
        
        for (ID n4 : edge_nbs){
            if (has_edge(n1, n2, n4) && lexic_comp_3(n1, n2, n3, n1, n2, n4))
                continue;
            if (has_edge(n1, n3, n4) && lexic_comp_3(n1, n2, n3, n1, n3, n4))
                continue;
            if (has_edge(n2, n3, n4) && lexic_comp_3(n1, n2, n3, n2, n3, n4))
                continue;

            auto bitcode {get_bitcode(n1, n2, n3, n4)};

            ++counts[motif_4_lookup[bitcode]];
            
            bitcode &= ~( (1u << 6) | (1u << 7) | (1u << 8) | (1u << 9));
            Mtype motif = {motif_4_lookup[bitcode]};
            if (motif != INV) { --counts[motif]; }
        }
    }
}

void Hypergraph::_motif_count_4_size_4(std::array<long long, 171>& counts) const
{
    for (const auto &edge : m_edges){
        if (edge.size() != 4) { continue; }
        
        unsigned int bitcode {get_bitcode(edge[0], edge[1], edge[2], edge[3])};

        Mtype motif = {motif_4_lookup[bitcode]};
        if (motif != INV)
            --counts[motif];
        
        bitcode |= (1u << 10);
        ++counts[motif_4_lookup[bitcode]];
    }
}


Hypergraph read_hyperedge_list(const std::string& filename)
{
	std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: could not open file: " + filename);
    }

	Hypergraph H {};
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;  // skip empty lines

        std::istringstream iss(line);
        std::vector<ID> edge;
        ID node;

        while (iss >> node) {
            edge.push_back(node);
        }

        if (!edge.empty() && !H.has_edge(edge)) {
            H.add_edge(edge);
        }
    }
	return H;
}
