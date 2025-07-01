#pragma once
#include <memory>
#include <optional>

template <typename Item, typename Priority, typename Key> class UpdatableQueue {
public:
  virtual ~UpdatableQueue() = default;

  virtual bool push(std::unique_ptr<Item> &&item_ptr,
                    const Priority &priority) = 0;

  virtual std::optional<std::unique_ptr<Item>> pop() = 0;

  virtual bool remove(const Key &key) = 0;

  virtual bool update(std::unique_ptr<Item> &&item_ptr,
                      const Priority &newPriority) = 0;

  virtual bool set(std::unique_ptr<Item> &&item_ptr,
                   const Priority &newPriority) = 0;

  virtual bool empty() const = 0;

  virtual size_t size() const = 0;
};
