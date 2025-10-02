#include "channel_util.hpp"

#include "global.hpp"
#include "io_channel_impl.hpp"

namespace ZebraGRAM::Channel {

void ChannelType::send_data(const void* data, uint64_t size) {
  switch (type) {
    case MEM_IO:
      static_cast<emp::MemIO*>(internal_channel)->send_data(data, size);
      break;
    case NET_IO:
      static_cast<emp::NetIO*>(internal_channel)->send_data(data, size);
      break;
    case FILE_IO:
      static_cast<emp::FileIO*>(internal_channel)->send_data(data, size);
      break;
    case HIGH_SPEED_NET_IO:
      static_cast<emp::HighSpeedNetIO*>(internal_channel)
          ->send_data(data, size);
      break;
    default:
      Assert(false);
  }
}

void ChannelType::recv_data(void* data, uint64_t size) {
  switch (type) {
    case MEM_IO:
      static_cast<emp::MemIO*>(internal_channel)->recv_data(data, size);
      break;
    case NET_IO:
      static_cast<emp::NetIO*>(internal_channel)->recv_data(data, size);
      break;
    case FILE_IO:
      static_cast<emp::FileIO*>(internal_channel)->recv_data(data, size);
      break;
    case HIGH_SPEED_NET_IO:
      static_cast<emp::HighSpeedNetIO*>(internal_channel)
          ->recv_data(data, size);
      break;
    default:
      Assert(false);
  }
}

bool is_channel_id(int id) {
  return ChannelManager::get_instance().is_channel_id(id);
}

ssize_t write(int fid, const void* data, ssize_t size) {
  ChannelManager& channel_manager = ChannelManager::get_instance();
  Assert(channel_manager.is_channel_id(fid));
  ChannelType channel = channel_manager.get_channel(fid);
  channel.send_data(data, size);
  return size;
}

int close(int fid) {
  ChannelManager& channel_manager = ChannelManager::get_instance();
  Assert(channel_manager.is_channel_id(fid));
  return channel_manager.delete_channel(fid);
}

}  // namespace ZebraGRAM::Channel
