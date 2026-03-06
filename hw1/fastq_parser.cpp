#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " <fastq format file path>\n";
    return 1;
  }

  std::ifstream fastq_file(argv[1]);
  if (!fastq_file.is_open()) {
    std::cerr << "could not open file " << argv[1] << '\n';
    return 1;
  }

  std::string line;
  uint64_t total_reads = 0;
  uint64_t total_len = 0;
  uint64_t min_len = std::numeric_limits<uint64_t>::max();
  uint64_t max_len = 0;
  uint64_t gc_counter = 0;

  const int PHRED_MAGIC_NUM = 33;
  uint64_t q_sum = 0;
  uint64_t long_reads_sum = 0;

  while (std::getline(fastq_file, line)) {
    std::string sequence;
    if (std::getline(fastq_file, sequence)) {
      uint64_t current_len = sequence.length();

      total_reads++;
      total_len += current_len;

      // min & max
      if (current_len < min_len)
        min_len = current_len;
      if (current_len > max_len)
        max_len = current_len;

      // gc-quality
      gc_counter += std::count_if(sequence.begin(), sequence.end(),
                                  [](char c) { return c == 'G' || c == 'C'; });

      // ignore comment
      fastq_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      // basecall qualities
      std::string basecall_qualities;
      std::getline(fastq_file, basecall_qualities);
      if (basecall_qualities.size() >= 10) {
        q_sum += (basecall_qualities[9] - PHRED_MAGIC_NUM);
        ++long_reads_sum;
      }
    }
  }

  fastq_file.close();

  if (total_reads == 0) {
    std::cout << "no stats for an empty file, sorry :(\n";
    return 0;
  }

  // final count of stats
  uint64_t average_len =
      std::round(static_cast<double>(total_len) / total_reads);
  uint64_t gc_content = (static_cast<double>(gc_counter) * 100) / total_len;
  uint64_t phred_quality =
      std::round(static_cast<double>(q_sum) / long_reads_sum);

  std::cout << "stats:\n";
  std::cout << "total reads: " << total_reads << '\n';
  std::cout << "minimal length of read: " << min_len << '\n';
  std::cout << "average length of read: " << average_len << '\n';
  std::cout << "maximal length of read: " << max_len << '\n';
  std::cout << "GC-content: " << gc_content << '\n';
  std::cout << "phred quality: " << phred_quality << '\n';

  return 0;
}
