#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <algorithm>

#include "utilities/unordered_dense.h"


using ID = unsigned int;


struct ArrayHash {
    template <typename T, size_t N>
    size_t operator()(const std::array<T, N>& a) const
	{
        size_t hash = 0;

        for (const auto& x : a) {
            hash ^= std::hash<T>{}(x) 
                  + 0x9e3779b9 
                  + (hash << 6) 
                  + (hash >> 2);
        }

        return hash;
    }
};


template <typename T>
bool lexic_comp_3(T n1, T n2, T n3, T m1, T m2, T m3)
{
	if (m2 < m1) { std::swap(m1,m2); }
	if (m3 < m2) { std::swap(m2,m3); }
	if (m2 < m1) { std::swap(m1,m2); }

	if (n1 != m1) { return n1 < m1; }
	if (n2 != m2) { return n2 < m2; }
	return n3 < m3;
}



class Hypergraph
{
private:
	ankerl::unordered_dense::map<ID, std::vector<ID>> m_nodes {};
	std::vector<std::vector<ID>> m_edges {};

	ankerl::unordered_dense::set<std::array<ID, 2>, ArrayHash> m_edge2_ids {};
    ankerl::unordered_dense::set<std::array<ID, 3>, ArrayHash> m_edge3_ids {};
    ankerl::unordered_dense::set<std::array<ID, 4>, ArrayHash> m_edge4_ids {};

	ID m_next_id {};

	void add_to_index(ID edge_id);

public:
	void add_node(ID n);

	void add_edge(std::vector<ID> nodes);

	const std::vector<ID>& get_edge(ID e) const { return m_edges.at(e); }

	auto node_count() const { return m_nodes.size(); }

	auto edge_count() const { return m_edges.size(); }

	inline bool has_edge(const std::vector<ID>& edge) const
    {
        if (edge.size() == 2) { return has_edge(edge[0], edge[1]); }
        if (edge.size() == 3) { return has_edge(edge[0], edge[1], edge[2]); }
        if (edge.size() == 4) { return has_edge(edge[0], edge[1], edge[2], edge[3]); }
        return false;
    }

	inline bool has_edge(ID n1, ID n2) const
	{
        if (n1 < n2) { return m_edge2_ids.contains({n1,n2}); }
        else { return m_edge2_ids.contains({n2,n1}); }
	}

    inline bool has_edge(ID n1, ID n2, ID n3) const
	{
        if (n1 > n2) { std::swap(n1, n2); }
        if (n2 > n3) { std::swap(n2, n3); }
        if (n1 > n2) { std::swap(n1, n2); }
        return m_edge3_ids.contains({n1,n2,n3});
	}

    inline bool has_edge(ID n1, ID n2, ID n3, ID n4) const
	{
        std::array<ID, 4> edge {n1,n2,n3,n4};
        std::sort(edge.begin(), edge.end());
        return m_edge4_ids.contains(edge);
	}

	std::vector<std::vector<ID>> edge_list() const { return m_edges; }

	int degree(ID n, size_t edge_size = 0) const;

	size_t size(ID e) const { return m_edges[e].size(); }

	auto neighbors(ID n, size_t edge_size = 0) const;

	std::vector<int> degree_list() const;

	std::vector<int> neighbor_count() const;

    std::vector<size_t> edge_size_count() const;


	std::array<long long, 6> motif_count_3() const;

	std::array<long long, 171> motif_count_4() const ;

	inline unsigned int get_bitcode(ID n1, ID n2, ID n3, ID n4) const
	{
		unsigned int bitcode {};
		if (has_edge(n1,n2)) { bitcode |= (1u << 0); }
		if (has_edge(n1,n3)) { bitcode |= (1u << 1); }
		if (has_edge(n1,n4)) { bitcode |= (1u << 2); }
		if (has_edge(n2,n3)) { bitcode |= (1u << 3); }
		if (has_edge(n2,n4)) { bitcode |= (1u << 4); }
		if (has_edge(n3,n4)) { bitcode |= (1u << 5); }
		if (has_edge(n1,n2,n3)) { bitcode |= (1u << 6); }
		if (has_edge(n1,n2,n4)) { bitcode |= (1u << 7); }
		if (has_edge(n1,n3,n4)) { bitcode |= (1u << 8); }
		if (has_edge(n2,n3,n4)) { bitcode |= (1u << 9); }

		return bitcode;
	}

	void _motif_count_4_size_2(std::array<long long, 171>& counts) const;

	void _motif_count_4_size_3(std::array<long long, 171>& counts) const;

	void _motif_count_4_size_4(std::array<long long, 171>& counts) const;
};


Hypergraph read_hyperedge_list(const std::string& filename);
