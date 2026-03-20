#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

struct Edge {
  std::string sequence;
  double coverage;
  std::string target_node;
};

class DeBruijnGraph {
private:
  int k;
  std::unordered_map<std::string, std::vector<Edge>> adj;
  std::unordered_map<std::string, int> in_degree;

public:
  DeBruijnGraph(int k_val) : k(k_val) {}

  void add_sequence(const std::string &seq) {
    if (seq.length() <= k)
      return;
    for (size_t i = 0; i <= seq.length() - (k + 1); ++i) {
      std::string u = seq.substr(i, k);
      std::string v = seq.substr(i + 1, k);
      std::string edge_seq = seq.substr(i, k + 1);

      bool found = false;
      for (auto &edge : adj[u]) {
        if (edge.target_node == v) {
          edge.coverage += 1.0;
          found = true;
          break;
        }
      }
      if (!found) {
        adj[u].push_back({edge_seq, 1.0, v});
        in_degree[v]++;
      }
    }
  }

  void clean_graph(double min_coverage) {
    for (auto &[u, edges] : adj) {
      edges.erase(remove_if(edges.begin(), edges.end(),
                            [&](const Edge &e) {
                              if (e.coverage < min_coverage) {
                                in_degree[e.target_node]--;
                                return true;
                              }
                              return false;
                            }),
                  edges.end());
    }
  }

  void compress() {
    bool simplified = true;
    while (simplified) {
      simplified = false;
      for (auto it = adj.begin(); it != adj.end(); ++it) {
        std::string u = it->first;
        if (adj[u].size() == 1) {
          std::string v = adj[u][0].target_node;
          if (u != v && in_degree[v] == 1 && adj.count(v) &&
              adj[v].size() == 1) {
            Edge &e1 = adj[u][0];
            Edge &e2 = adj[v][0];

            e1.sequence += e2.sequence.substr(k);
            e1.coverage = (e1.coverage + e2.coverage) / 2.0;
            e1.target_node = e2.target_node;

            adj.erase(v);
            in_degree.erase(v);
            simplified = true;
            break;
          }
        }
      }
    }
  }

  void save_gfa(const std::string &filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
      std::cerr << "error: could not create .gfa file\n";
      return;
    }
    out << "H\tVN:Z:1.0\n";

    // segment info
    struct Seg {
      int id;
      std::string seq;
      double cov;
      std::string u, v; // s and t node
    };

    std::vector<Seg> segments;
    std::unordered_map<std::string, std::vector<int>> incoming;
    std::unordered_map<std::string, std::vector<int>> outgoing;

    int seg_id = 1;
    for (auto const &[u, edges] : adj) {
      for (auto const &e : edges) {
        int id = seg_id++;
        segments.push_back({id, e.sequence, e.coverage, u, e.target_node});
        outgoing[u].push_back(id);
        incoming[e.target_node].push_back(id);
      }
    }

    for (auto const &s : segments) {
      out << "S\t" << s.id << "\t" << s.seq << "\tKC:i:" << (int)s.cov << "\n";
    }

    for (auto const &[node, inc_ids] : incoming) {
      if (outgoing.count(node)) {
        for (int in_id : inc_ids) {
          for (int out_id : outgoing[node]) {
            out << "L\t" << in_id << "\t+\t" << out_id << "\t+\t" << k << "M\n";
          }
        }
      }
    }
    out.close();
  }

  void save_fasta(const std::string &filename) {
    std::ofstream out(filename);
    if (!out.is_open())
      return;
    int count = 0;
    for (auto const &[u, edges] : adj) {
      for (auto const &e : edges) {
        out << ">contig_" << ++count << "_len_" << e.sequence.length()
            << "_cov_" << e.coverage << "\n";
        out << e.sequence << "\n";
      }
    }
    out.close();
  }
};

std::vector<std::string> parse_data(const std::string &path) {
  std::vector<std::string> reads;
  std::ifstream f(path);
  if (!f.is_open()) {
    std::cerr << "error: file " << path << " not found!";
    return reads;
  }
  std::string line, seq;
  bool is_fastq = (path.find(".fastq") != std::string::npos ||
                   path.find(".fq") != std::string::npos);

  while (getline(f, line)) {
    if (line.empty())
      continue;
    if (line[0] == '>' || line[0] == '@') {
      if (getline(f, seq)) {
        reads.push_back(seq);
        if (is_fastq) {
          getline(f, line);
          getline(f, line);
        }
      }
    }
  }
  return reads;
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout << "usage: " << argv[0] << " <input> <k> <prefix> [min_cov]\n";
    return 1;
  }

  std::string input = argv[1];
  int k = std::stoi(argv[2]);
  std::string prefix = argv[3];
  double min_cov = (argc > 4) ? std::stod(argv[4]) : 1.0;

  std::cout << "starting..." << std::endl;
  auto reads = parse_data(input);
  if (reads.empty())
    return 1;

  std::cout << "building graph...\n";
  DeBruijnGraph dbg(k);
  for (const auto &r : reads)
    dbg.add_sequence(r);

  if (min_cov > 1.0) {
    std::cout << "cleaning (cov < " << min_cov << ")...\n";
    dbg.clean_graph(min_cov);
  }

  std::cout << "compressing...\n";
  dbg.compress();

  std::cout << "saving...\n";
  dbg.save_gfa(prefix + ".gfa");
  dbg.save_fasta(prefix + ".fasta");

  std::cout << "done!\n";
  std::cout << "results in files: " << std::filesystem::current_path() << "/"
            << prefix << ".gfa\n";

  return 0;
}
