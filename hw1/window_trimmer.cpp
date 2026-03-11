#include "window_trimmer.h"

// returns first index of trimmed part
uint64_t window_trimmer::find_trimmed_index(const std::string_view &qual_string,
                                            uint16_t width,
                                            double avg_quality) {
  uint8_t baseline = width * avg_quality;

  // string shorter than window -> invalid
  if (qual_string.size() < width) {
    return 0;
  }

  uint32_t wind_qual = 0;

  // count first window, so for the next windows
  // it is only needed to subtract tail and add head
  for (int i = 0; i != width; ++i) {
    wind_qual += (qual_string[i] - PHRED_ASCII_BASE);
  }

  // if first window is lower or equals than baseline,
  // then this read is trimmed
  if (wind_qual <= baseline) {
    return 0;
  }

  // main cycle
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
  // searching for the "good" prefix
  if (broken) {
    window_end = i + width - 1;
    while (window_end > i - 1 &&
           qual_string[window_end] < avg_quality + PHRED_ASCII_BASE) {
      --window_end;
    }
  }
  return window_end + 1;
}