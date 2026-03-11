#include "window_trimmer.h"

uint64_t window_trimmer::find_trimmed_index(const std::string_view &qual_string,
                                            int width) {
  uint8_t baseline = width * 30;
  if (qual_string.size() < width) {
    return 0;
  }

  uint32_t wind_qual = 0;

  for (int i = 0; i != width; ++i) {
    wind_qual += (qual_string[i] - PHRED_ASCII_BASE);
  }

  if (wind_qual <= baseline) {
    return 0;
  }

  uint32_t i = 1;
  bool broken = false;
  for (; i < qual_string.size() - width + 1; ++i) {
    wind_qual = wind_qual - qual_string[i - 1] + qual_string[i + width - 1];
    if (wind_qual <= baseline) {
      broken = true;
      break;
    }
  }
  uint32_t window_end = qual_string.size() - 1;
  if (broken) {
    window_end = i + width - 1;
    while (window_end > i - 1 && qual_string[window_end] < 63) {
      --window_end;
    }
  }
  return window_end + 1;
}