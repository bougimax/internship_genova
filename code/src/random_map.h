#include <iostream>
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <stdexcept>

#pragma once

template<typename K, typename V>
class RandomMap {
private:
    std::unordered_map<K, V> data;          
    std::vector<K> keys;                    
    std::unordered_map<K, size_t> indices; 

public:
    RandomMap() {
        std::srand(std::time(nullptr));
    }

    void insert(const K& key, const V& value) {
        if (data.find(key) == data.end()) {
            keys.push_back(key);
            indices[key] = keys.size() - 1;
        }
        data[key] = value;
    }

    std::pair<K, V> random_pop() {
        if (data.empty()) {
            throw std::out_of_range("pop from empty map");
        }

        int idx = std::rand() % keys.size();
        K key = keys[idx];

        // Swap dans keys
        K last_key = keys.back();
        std::swap(keys[idx], keys.back());

        // Update indices après swap
        indices[last_key] = idx;
        indices.erase(key);

        keys.pop_back();

        V value = data[key];
        data.erase(key);

        return {key, value};
    }

    void erase(const K& key) {
        auto it = data.find(key);
        if (it == data.end()) {
            return;
        }

        size_t idx = indices[key];
        K last_key = keys.back();

        // Swap la clé à supprimer avec la dernière dans keys
        std::swap(keys[idx], keys.back());

        // Mettre à jour l’index du dernier élément
        indices[last_key] = idx;

        // Supprimer le dernier élément
        keys.pop_back();
        indices.erase(key);
        data.erase(key);
    }

    const V& get(const K& key) const {
        auto it = data.find(key);
        if (it == data.end()) {
            throw std::invalid_argument("key not found");
        }
        return it->second;
    }

    bool contains(const K& key) const {
        return data.find(key) != data.end();
    }

    void print() const {
        std::cout << "{ ";
        for (const auto& k : keys) {
            std::cout << k << ": " << data.at(k) << ", ";
        }
        std::cout << "}" << std::endl;
    }

    bool empty() const {
        return data.empty();
    }

    size_t size() const {
        return data.size();
    }
};

