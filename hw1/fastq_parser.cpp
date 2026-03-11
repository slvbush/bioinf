#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "window_trimmer.h"

int main(int argc, char *argv[]) {
  // usage check
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " <fastq format file path>\n";
    return 1;
  }

  // open file, check
  std::ifstream fastq_file(argv[1]);
  if (!fastq_file.is_open()) {
    std::cerr << "could not open file " << argv[1] << '\n';
    return 1;
  }

  // all the variables for statistics
  uint64_t total_reads = 0;
  uint64_t total_len = 0;
  uint64_t min_len = std::numeric_limits<uint64_t>::max();
  uint64_t max_len = 0;
  uint64_t gc_counter = 0;

  const uint8_t PHRED_ASCII_BASE = 33;
  const uint8_t PHRED_IDX = 9;
  uint64_t qual_sum = 0;
  uint64_t long_reads_ctr = 0;

  uint64_t trimmed_counter = 0;

  uint64_t total_reads_left = 0;
  uint64_t total_len_left = 0;
  uint64_t min_len_left = std::numeric_limits<uint64_t>::max();
  uint64_t max_len_left = 0;

  uint64_t long_counter = 0;

  // read 1st string - name
  std::string line;
  while (std::getline(fastq_file, line)) {
    // read 2nd string - read
    std::string sequence;
    if (std::getline(fastq_file, sequence)) {

      // avg
      ++total_reads;
      uint64_t curr_len = sequence.length();
      total_len += curr_len;

      // min & max
      min_len = std::min(min_len, curr_len);
      max_len = std::max(max_len, curr_len);

      // gc-quality
      gc_counter += std::count_if(sequence.begin(), sequence.end(),
                                  [](char c) { return c == 'G' || c == 'C'; });

      // ignore 3d string - comment
      fastq_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      // read 4th string - basecall qualities
      std::string basecall_qualities;
      std::getline(fastq_file, basecall_qualities);

      // phred quality on position 10
      if (basecall_qualities.size() > PHRED_IDX) {
        qual_sum += (basecall_qualities[PHRED_IDX] - PHRED_ASCII_BASE);
        ++long_reads_ctr;
      }

      // window trimming + stats & length trimming after
      uint64_t tr_idx = window_trimmer::find_trimmed_index(
          basecall_qualities.substr(), 5, 30);
      if (tr_idx == 0) {
        ++trimmed_counter;
      } else {
        min_len_left = std::min(min_len_left, tr_idx);
        max_len_left = std::max(max_len_left, tr_idx);
        total_reads_left++;
        total_len_left += tr_idx;
        if (tr_idx >= 60) {
          ++long_counter;
        }
      }
    }
  }

  if (total_reads == 0) {
    std::cout << "no stats for an empty file, sorry :(\n";
    return 0;
  }

  // final count of stats
  uint64_t avg_len = std::round(static_cast<double>(total_len) / total_reads);
  double gc_content = static_cast<double>(gc_counter * 100) / total_len;
  uint64_t phred_quality =
      std::round(static_cast<double>(qual_sum) / long_reads_ctr);
  uint64_t avg_len_after_trimming =
      std::round(static_cast<double>(total_len_left) / total_reads_left);

  // outputs
  std::cout << "STATS:\n";

  std::cout << "total reads: " << total_reads << '\n';
  std::cout << "minimal length of read: " << min_len << '\n';
  std::cout << "maximal length of read: " << max_len << '\n';
  std::cout << "average length of read: " << avg_len << '\n';

  std::cout << '\n';

  std::cout << "GC-content: " << gc_content << '\n';
  std::cout << "phred quality: " << phred_quality << '\n';

  std::cout << '\n';

  std::cout << "trimmed with window: " << trimmed_counter << '\n';
  std::cout << "minimal length of read left after trimming: " << min_len_left
            << '\n';
  std::cout << "maximal length of read left after trimming: " << max_len_left
            << '\n';
  std::cout << "average length of read left after trimming: "
            << avg_len_after_trimming << '\n';

  std::cout << '\n';

  std::cout << "amount of reads left after second trimming with length 60: "
            << long_counter << '\n';
  return 0;
}
