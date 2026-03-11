#pragma once

#include <cstdint>
#include <string>

class window_trimmer {
public:
  static uint64_t find_trimmed_index(const std::string_view &qual, int width);

private:
  static constexpr uint8_t PHRED_ASCII_BASE = 33;
};
