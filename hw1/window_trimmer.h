#pragma once

#include <cstdint>
#include <string>

class window_trimmer {
public:
  // returns first index of trimmed part
  static uint64_t find_trimmed_index(const std::string_view &qual,
                                     uint16_t width, double avg_quality);

private:
  static constexpr uint8_t PHRED_ASCII_BASE = 33;
};
