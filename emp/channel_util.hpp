#pragma once
#include <sys/resource.h>

#include <vector>

#include "../global.hpp"
#include "../utils/util.hpp"
namespace PicoGRAM::Channel {

/**
 * @brief A wrapper for different types of IO channels
 *
 */
struct ChannelType {
 private:
  void* internal_channel = nullptr;
  IOType type = INVALID_IO;

 public:
  ChannelType() = default;

  ChannelType(void* channel, IOType type)
      : internal_channel(channel), type(type) {}

  /**
   * @brief Send data through the channel
   *
   * @param data the data to send
   * @param size the size of the data in bytes
   */
  void send_data(const void* data, uint64_t size);

  /**
   * @brief Receive data through the channel
   *
   * @param data the buffer to store the received data
   * @param size the size of the data in bytes
   */
  void recv_data(void* data, uint64_t size);

  /**
   * @brief Check if the channel is set
   *
   * @return true if the channel is set
   */
  bool is_set() const { return internal_channel != nullptr; }

  /**
   * @brief Unset the channel. Does not free the resource of the channel.
   *
   */
  void unset() {
    internal_channel = nullptr;
    type = INVALID_IO;
  }
  bool operator==(const ChannelType& other) const {
    return internal_channel == other.internal_channel && type == other.type;
  }
  bool operator!=(const ChannelType& other) const { return !(*this == other); }
};

/**
 * @brief A singleton class to manage channels, so that the client only needs to
 * keep track of the channel id, similar to a file descriptor.
 *
 */
struct ChannelManager {
 private:
  int start_id;  // the starting id for the channels, which is the maximum
                 // number of file descriptors allowed by the OS
  std::vector<ChannelType> channels;

  ChannelManager() {
    struct rlimit limit;
    getrlimit(RLIMIT_NOFILE, &limit);
    start_id = limit.rlim_max;
  }

 public:
  /**
   * @brief Get the singleton instance of the ChannelManager
   *
   * @return ChannelManager&
   */
  static ChannelManager& get_instance() {
    static ChannelManager instance;
    return instance;
  }

  /**
   * @brief Check if the id is a channel id
   *
   * @param id the id to check
   * @return true if the id is a channel id
   * @return false if the id is a file descriptor
   */
  bool is_channel_id(int id) const { return id >= start_id; }

  /**
   * @brief Add a channel to the manager
   *
   * @param channel the channel to add
   * @return int the id of the channel
   */
  int add_channel(ChannelType channel) {
    Assert(channels.size() < 1UL << 20);
    // first check if channel already exists
    for (int idx = 0; idx < (int)channels.size(); ++idx) {
      if (channels[idx] == channel) {
        return idx + start_id;
      }
    }
    int idx = 0;
    for (; idx < (int)channels.size(); ++idx) {
      if (!channels[idx].is_set()) {
        channels[idx] = channel;
        return idx + start_id;
      }
    }
    channels.push_back(channel);
    return idx + start_id;
  }

  /**
   * @brief Get the channel by id
   *
   * @param id the id of the channel
   * @return ChannelType the channel
   */
  ChannelType get_channel(int id) {
    Assert(is_channel_id(id));
    int idx = id - start_id;
    Assert((uint64_t)idx < channels.size());
    Assert(channels[idx].is_set());
    return channels[idx];
  }

  /**
   * @brief Delete the channel by id. Does not free the resource of the channel.
   *
   * @param id the id of the channel
   * @return int 0 if successful, -1 if the channel does not exist
   */
  int delete_channel(int id) {
    if (!is_channel_id(id)) {
      return -1;
    }
    int idx = id - start_id;
    if (idx >= (int)channels.size()) {
      return -1;
    }
    if (!channels[idx].is_set()) {
      return -1;
    }

    channels[idx].unset();
    return 0;
  }
};

bool is_channel_id(int id);

/**
 * @brief Write data to a file descriptor or send it through a channel,
 * depending on the fid
 *
 * @param fid the file descriptor or channel id
 * @param data the data to write
 * @param size the size of the data in bytes
 * @return ssize_t the number of bytes written
 */
ssize_t write(int fid, const void* data, ssize_t size);

/**
 * @brief Close a file descriptor or delete a channel, depending on the fid
 *
 * @param fid the file descriptor or channel id
 * @return int 0 if successful, -1 if the file descriptor or channel does not
 * exist
 */
int close(int fid);
}  // namespace PicoGRAM::Channel
