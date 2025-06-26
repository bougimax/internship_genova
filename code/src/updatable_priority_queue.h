#pragma once
#include "updatable_queue_template.h"
#include <memory>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

template <typename Item, typename Priority, typename Key>
struct priority_queue_node {
  Priority priority;
  std::unique_ptr<Item> data;
  Key key;
  priority_queue_node(std::unique_ptr<Item> &&item_ptr,
                      const Priority &priority)
      : priority(priority), data(std::move(item_ptr)), key(data->id) {}
  friend bool operator<(const priority_queue_node &pqn1,
                        const priority_queue_node &pqn2) {
    return pqn1.priority < pqn2.priority;
  }
  friend bool operator>(const priority_queue_node &pqn1,
                        const priority_queue_node &pqn2) {
    return pqn1.priority > pqn2.priority;
  }
};

/** Key has to be an uint value (convertible to size_t)
 * This is a max heap (max is on top), to match stl's pQ */
template <typename Item, typename Priority, typename Key>
class UpdatablePriorityQueue : public UpdatableQueue<Item, Priority, Key> {
protected:
  std::unordered_map<Key, size_t> id_to_heappos;
  std::vector<priority_queue_node<Item, Priority, Key>> heap;

public:
  UpdatablePriorityQueue() {}

  bool empty() const override { return heap.empty(); }
  std::size_t size() const override { return heap.size(); }

  std::optional<std::unique_ptr<Item>> pop() override {
    if (size() == 0)
      return std::nullopt;
    priority_queue_node<Item, Priority, Key> ret = std::move(*heap.begin());
    id_to_heappos[ret.key] = -1;
    if (size() > 1) {
      *heap.begin() = std::move(*(heap.end() - 1));
      id_to_heappos[heap.front().key] = 0;
    }
    heap.pop_back();
    sift_down(0);
    if (ret.priority < 0)
      return std::nullopt;
    return std::make_optional(std::move(ret.data));
  }

  /** Sets the priority for the given key. If not present, it will be added,
   * otherwise it will be updated Returns true if the priority was changed.
   * */
  bool set(std::unique_ptr<Item> &&item_ptr,
           const Priority &priority) override {
    if (id_to_heappos.contains(item_ptr->id) &&
        id_to_heappos[item_ptr->id] <
            ((size_t)-2)) // This key is already in the pQ
      return update(std::move(item_ptr), priority);
    else
      return push(std::move(item_ptr), priority);
  }

  bool remove(const Key &id) override {
    if (!id_to_heappos.contains(id))
      return false;
    size_t heappos = id_to_heappos[id];
    if (heappos >= ((size_t)-2))
      return false;
    Priority &priority = heap[heappos].priority;
    priority = -1;
    sift_down(heappos);
    return true;
  }

  /** Returns true if the key was not inside and was added, otherwise does
   * nothing and returns false If the key was remembered and only_if_unknown is
   * true, does nothing and returns false
   * */
  bool push(std::unique_ptr<Item> &&item_ptr,
            const Priority &priority) override {
    if (id_to_heappos.contains(item_ptr->id))
      return false;
    // otherwise we have id_to_heappos[key] = -1, unseen key
    size_t n = heap.size();
    id_to_heappos[item_ptr->id] =
        n; // For consistency in the case where nothing moves (early return)
    heap.emplace_back(std::move(item_ptr), priority);
    sift_up(n);
    return true;
  }

  /** Returns true if the key was already inside and was updated, otherwise does
   * nothing and returns false */
  bool update(std::unique_ptr<Item> &&item_ptr,
              const Priority &new_priority) override {
    if (!id_to_heappos.contains(item_ptr->id))
      return false;
    size_t heappos = id_to_heappos[item_ptr->id];
    if (heappos >= ((size_t)-2))
      return false;
    Priority &priority = heap[heappos].priority;
    std::unique_ptr<Item> &data = heap[heappos].data;
    data = std::move(item_ptr);
    if (new_priority > priority) {
      priority = new_priority;
      sift_up(heappos);
      return true;
    } else if (new_priority < priority) {
      priority = new_priority;
      sift_down(heappos);
      return true;
    }
    return false;
  }

private:
  void extend_ids(Key k) {
    size_t new_size = k + 1;
    if (id_to_heappos.size() < new_size)
      id_to_heappos.resize(new_size, -1);
  }

  void sift_down(size_t heappos) {
    size_t len = heap.size();
    size_t child = heappos * 2 + 1;
    if (len < 2 || child >= len)
      return;
    if (child + 1 < len && heap[child + 1] > heap[child])
      ++child; // Check whether second child is higher
    if (!(heap[child] > heap[heappos]))
      return; // Already in heap order

    priority_queue_node<Item, Priority, Key> val = std::move(heap[heappos]);
    do {
      heap[heappos] = std::move(heap[child]);
      id_to_heappos[heap[heappos].key] = heappos;
      heappos = child;
      child = 2 * child + 1;
      if (child >= len)
        break;
      if (child + 1 < len && heap[child + 1] > heap[child])
        ++child;
    } while (heap[child] > val);
    heap[heappos] = std::move(val);
    id_to_heappos[heap[heappos].key] = heappos;
  }

  void sift_up(size_t heappos) {
    size_t len = heap.size();
    if (len < 2 || heappos <= 0)
      return;
    size_t parent = (heappos - 1) / 2;
    if (!(heap[heappos] > heap[parent]))
      return;
    priority_queue_node<Item, Priority, Key> val = std::move(heap[heappos]);
    do {
      heap[heappos] = std::move(heap[parent]);
      id_to_heappos[heap[heappos].key] = heappos;
      heappos = parent;
      if (heappos <= 0)
        break;
      parent = (parent - 1) / 2;
    } while (val > heap[parent]);
    heap[heappos] = std::move(val);
    id_to_heappos[heap[heappos].key] = heappos;
  }
};
