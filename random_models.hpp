#pragma once

#include <vector>

#include "hypergraph_motif.hpp"


struct VectorHash {
	template <typename T>
	size_t operator()(const std::vector<T> &v) const
	{
		size_t hash = v.size();
		for (auto &i : v) {
			hash ^= std::hash<ID>{}(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		}
		return hash;
	}
};

Hypergraph sample_sizepreserving(const Hypergraph& H_orig);

Hypergraph sample_configuration(const Hypergraph& H_orig, int swaps);

std::vector<double> motif_abundance(const Hypergraph& H, int motif_size, int model=1, int samples=100, int swaps=-1, double epsilon=4.0);
// model 0: size-preserving
// model 1: configuration

