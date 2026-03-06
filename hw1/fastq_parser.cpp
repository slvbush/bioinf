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

  while (std::getline(fastq_file, line)) {
    std::string sequence;
    if (std::getline(fastq_file, sequence)) {
      uint64_t current_len = sequence.length();

      total_reads++;
      total_len += current_len;

      if (current_len < min_len)
        min_len = current_len;
      if (current_len > max_len)
        max_len = current_len;

      gc_counter += std::count_if(sequence.begin(), sequence.end(),
                                  [](char c) { return c == 'G' || c == 'C'; });

      fastq_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      fastq_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }

  fastq_file.close();

  if (total_reads == 0) {
    std::cout << "no stats for an empty file, sorry :(\n";
    return 0;
  }

  long long average_len =
      std::round(static_cast<double>(total_len) / total_reads);

  std::cout << "stats:\n";
  std::cout << "total reads: " << total_reads << '\n';
  std::cout << "minimal length of read: " << min_len << '\n';
  std::cout << "average length of read: " << average_len << '\n';
  std::cout << "maximal length of read: " << max_len << '\n';
  std::cout << "GC-content: "
            << (static_cast<double>(gc_counter) * 100) / total_len << '\n';

  return 0;
}
