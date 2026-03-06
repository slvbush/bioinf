#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " <fastq format file path>"
              << std::endl;
    return 1;
  }

  std::ifstream fastq_file(argv[1]);
  if (!fastq_file.is_open()) {
    std::cerr << "could not open file " << argv[1] << std::endl;
    return 1;
  }

  std::string line;
  long long total_reads = 0;
  long long total_length = 0;
  size_t min_length = std::numeric_limits<size_t>::max();
  size_t max_length = 0;

  while (std::getline(fastq_file, line)) {
    std::string sequence;
    if (std::getline(fastq_file, sequence)) {
      size_t current_len = sequence.length();

      total_reads++;
      total_length += current_len;

      if (current_len < min_length)
        min_length = current_len;
      if (current_len > max_length)
        max_length = current_len;

      std::string plus_line;
      std::string quality;
      std::getline(fastq_file, plus_line);
      std::getline(fastq_file, quality);
    }
  }

  fastq_file.close();

  if (total_reads == 0) {
    std::cout << "no stats for an empty file, sorry :(" << std::endl;
    return 0;
  }

  long long average_length =
      std::round(static_cast<double>(total_length) / total_reads);

  std::cout << "stats:\n";
  std::cout << "total reads: " << total_reads << '\n';
  std::cout << "minimal length of read: " << min_length << '\n';
  std::cout << "average length of read: " << average_length << '\n';
  std::cout << "maximal length of read: " << max_length << '\n';

  return 0;
}
