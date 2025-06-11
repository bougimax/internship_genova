#include "updatable_queue_template.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#pragma once

template <typename Item, typename Key>
class UpdatableRandomQueue : public UpdatableQueue<Item, double, Key> {
private:
  std::unordered_map<Key, std::unique_ptr<Item>> data;
  std::vector<Key> keys;
  std::unordered_map<Key, size_t> indices;

public:
  UpdatableRandomQueue() { std::srand(std::time(nullptr)); }

  bool push(std::unique_ptr<Item> &&item_ptr, const double &priority) override {
    if (data.find(item_ptr->id) == data.end()) {
      keys.push_back(item_ptr->id);
      indices[item_ptr->id] = keys.size() - 1;
    }
    data[item_ptr->id] = std::move(item_ptr);
    return true;
  }

  std::optional<std::unique_ptr<Item>> pop() override {
    if (data.empty()) {
      return std::nullopt;
    }

    int idx = std::rand() % keys.size();
    Key key = keys[idx];

    // Swap dans keys
    Key last_key = keys.back();
    std::swap(keys[idx], keys.back());

    // Update indices après swap
    indices[last_key] = idx;
    indices.erase(key);

    keys.pop_back();

    std::unique_ptr<Item> value = std::move(data[key]);
    data.erase(key);

    return std::make_optional(std::move(value));
  }

  bool remove(const Key &key) override {
    auto it = data.find(key);
    if (it == data.end()) {
      return false;
    }

    size_t idx = indices[key];
    Key last_key = keys.back();

    // Swap la clé à supprimer avec la dernière dans keys
    std::swap(keys[idx], keys.back());

    // Mettre à jour l’index du dernier élément
    indices[last_key] = idx;

    // Supprimer le dernier élément
    keys.pop_back();
    indices.erase(key);
    data.erase(key);

    return true;
  }

  bool update(std::unique_ptr<Item> &&item_ptr,
              const double &newPriority) override {
    if (data.find(item_ptr->id) != data.end()) {
      push(std::move(item_ptr), newPriority);
      return true;
    }
    return false;
  }

  bool set(std::unique_ptr<Item> &&item_ptr,
           const double &newPriority) override {
    return push(std::move(item_ptr), newPriority);
  }

  bool empty() const override { return data.empty(); }

  size_t size() const override { return keys.size(); }

private:
  bool contains(const Item &item) const {
    return data.find(item->id) != data.end();
  }
};
