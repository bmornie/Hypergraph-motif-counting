#pragma once

#include <random>
#include <chrono>
#include <vector>
#include <unordered_set>
#include <algorithm>


inline std::mt19937 generate()
{
    std::random_device rd {};

    // Create seed_seq with clock and 7 random numbers from std::random_device
    std::seed_seq ss{
        static_cast<std::seed_seq::result_type>(std::chrono::steady_clock::now().time_since_epoch().count()),
            rd(), rd(), rd(), rd(), rd(), rd(), rd() };

    return std::mt19937{ ss };
}

inline std::mt19937& rng()
{
    thread_local std::mt19937 gen {generate()};
    return gen;
}


template <typename T>
T rand_int(T lower, T upper)
{
    std::uniform_int_distribution<T> dist(lower, upper);
    return dist(rng());
}

template<typename T>
std::vector<T> rand_int(T lower, T upper, size_t k)
{
    std::unordered_set<T> result;
    std::uniform_int_distribution<T> dist(lower, upper);
    while (result.size() < k) { result.insert(dist(rng())); }
    return std::vector<T>(result.begin(), result.end());
}


template <typename T>
T pop_random(std::vector<T>& vec)
{
    std::uniform_int_distribution<size_t> dist(0, vec.size() - 1);
    size_t index = dist(rng());

    T value = vec[index];
    std::swap(vec[index], vec.back());
    vec.pop_back();

    return value;
}