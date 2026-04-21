#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>

#include "SBWT.hh"
#include "SeqIO/SeqIO.hh"
#include "SeqIO/buffered_streams.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetWT.hh"
#include "commands.hh"
#include "cxxopts.hpp"
#include "globals.hh"
#include "stdlib_printing.hh"
#include "variants.hh"

using namespace std;

using namespace sbwt;

// Assumes values of v are -1 or larger
template <typename writer_t>
inline void print_vector(const vector<int64_t>& v, writer_t& out) {
  // Fast manual integer-to-string conversion
  char buffer[32];
  char newline = '\n';
  for (int64_t x : v) {
    int64_t i = 0;
    if (x == -1) {
      buffer[0] = '1';
      buffer[1] = '-';
      i = 2;
    } else {
      while (x > 0) {
        buffer[i++] = '0' + (x % 10);
        x /= 10;
      }
    }
    std::reverse(buffer, buffer + i);
    buffer[i] = ' ';
    out.write(buffer, i + 1);
  }
  out.write(&newline, 1);
}

template <typename T>
std::vector<T> read_file(const char* filename, size_t offset = 0) {
  std::ifstream ifs(filename, std::ios::binary);

  const auto begin = ifs.tellg();
  ifs.seekg(0, std::ios::end);
  const auto end = ifs.tellg();
  const std::size_t len = (end - begin) / sizeof(T);
  ifs.seekg(0);

  std::vector<T> v(len, 0);

  for (std::size_t i = 0; i < len; ++i) {
    ifs.read(reinterpret_cast<char*>(v.data() + i), sizeof(T));
  }

  ifs.close();

  return v;
}

template <typename sbwt_t, typename writer_t>
int64_t run_file(const string& infile_p, const string& infile_c,
                 const string& outfile, const sbwt_t& sbwt) {
  writer_t writer(outfile);
  std::vector<char> query_chars = read_file<char>(infile_c.c_str());
  std::vector<int64_t> query_pos = read_file<int64_t>(infile_p.c_str());
  if (query_chars.size() != query_pos.size()) {
    std::cerr << "Character query file size: " << query_chars.size()
              << std::endl;
    std::cerr << "Position query file size: " << query_pos.size() << std::endl;
    throw std::runtime_error("Index and character files do not match in size");
  }
  auto s_rank_struct = sbwt.get_subset_rank_structure();
  auto n_nodes = sbwt.number_of_subsets();
  std::cerr << "Running " << query_chars.size()
            << " s-rank queries on SBWT with " << n_nodes << " nodes."
            << std::endl;
  int64_t total_micros = 0;
  int64_t number_of_queries = 0;

  uint64_t total_ranks = 0;
  size_t n = query_chars.size();
  vector<int64_t> out_buffer;
  for (size_t i = 0; i < n; ++i) {
    if (query_pos[i] < 0) {
      throw std::runtime_error("Negative query length found");
    }
    if (query_pos[i] >= n_nodes) {
      throw std::runtime_error(
          "Query position out of bounds: " + std::to_string(query_pos[i]) +
          " >= " + std::to_string(n_nodes));
    }
    total_ranks %= 1000000007;

    char c = query_chars[i];
    int64_t l = query_pos[i];

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, c);
    total_micros += cur_time_micros() - t0;
    number_of_queries++;
    total_ranks += ans;
    out_buffer.push_back(ans);
  }
  std::cerr << "Total ranks: " << total_ranks << std::endl;
  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);
  print_vector(out_buffer, writer);
  out_buffer.clear();

  std::cerr << "Running " << query_chars.size()
            << " s-rank queries with pair: for every srr query run additional "
               "query for the position close by (+5)."
            << std::endl;
  total_micros = 0;
  number_of_queries = 0;
  total_ranks = 0;
  n = query_chars.size();

  for (size_t i = 0; i < n; ++i) {
    if (query_pos[i] < 0) {
      throw std::runtime_error("Negative query length found");
    }
    if (query_pos[i] >= n_nodes) {
      throw std::runtime_error(
          "Query position out of bounds: " + std::to_string(query_pos[i]) +
          " >= " + std::to_string(n_nodes));
    }
    total_ranks %= 1000000007;

    char c = query_chars[i];
    int64_t l = query_pos[i];
    int64_t l2 = l + 5;
    if (l2 >= n_nodes) l2 = n_nodes - 1;

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, c);
    auto ans_2 = s_rank_struct.rank(l2, c);
    total_micros += cur_time_micros() - t0;
    number_of_queries += 2;
    total_ranks += ans + ans_2;
    out_buffer.push_back(ans);
  }
  std::cerr << "Total ranks: " << total_ranks << std::endl;
  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);
  print_vector(out_buffer, writer);
  out_buffer.clear();

  std::cerr << "Running " << query_chars.size()
            << " s-rank queries for all characters: for every srr query run "
               "additional run the position for all 4 characters A, C, G, T."
            << std::endl;
  total_micros = 0;
  number_of_queries = 0;
  total_ranks = 0;
  n = query_chars.size();

  for (size_t i = 0; i < n; ++i) {
    if (query_pos[i] < 0) {
      throw std::runtime_error("Negative query length found");
    }
    if (query_pos[i] >= n_nodes) {
      throw std::runtime_error(
          "Query position out of bounds: " + std::to_string(query_pos[i]) +
          " >= " + std::to_string(n_nodes));
    }
    total_ranks %= 1000000007;

    char c = query_chars[i];
    int64_t l = query_pos[i];
    int64_t l2 = l + 5;
    if (l2 >= n_nodes) l2 = n_nodes - 1;

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, 'A');
    auto ans_2 = s_rank_struct.rank(l, 'C');
    auto ans_3 = s_rank_struct.rank(l, 'G');
    auto ans_4 = s_rank_struct.rank(l, 'T');
    total_micros += cur_time_micros() - t0;
    number_of_queries += 4;
    total_ranks += ans + ans_2 + ans_3 + ans_4;
    out_buffer.push_back(ans);
  }
  std::cerr << "Total ranks: " << total_ranks << std::endl;
  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);
  print_vector(out_buffer, writer);
  out_buffer.clear();

  write_log("Now running on an array of random positions for character A",
            LogLevel::MAJOR);
  total_micros = 0;
  for (size_t i = 0; i < n; ++i) {
    if (query_pos[i] < 0) {
      throw std::runtime_error("Negative query length found");
    }
    if (query_pos[i] >= n_nodes) {
      throw std::runtime_error(
          "Query position out of bounds: " + std::to_string(query_pos[i]) +
          " >= " + std::to_string(n_nodes));
    }
    total_ranks %= 1000000007;

    char c = 'A';
    int64_t l = query_pos[i];

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, c);
    total_micros += cur_time_micros() - t0;
    total_ranks += ans;
  }
  std::cerr << "Total ranks: " << total_ranks << std::endl;
  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);

  write_log("Now running on an array of random positions for character C",
            LogLevel::MAJOR);

  total_ranks = 0;
  total_micros = 0;
  for (size_t i = 0; i < n; ++i) {
    if (query_pos[i] < 0) {
      throw std::runtime_error("Negative query length found");
    }
    if (query_pos[i] >= n_nodes) {
      throw std::runtime_error(
          "Query position out of bounds: " + std::to_string(query_pos[i]) +
          " >= " + std::to_string(n_nodes));
    }
    total_ranks %= 1000000007;

    char c = 'C';
    int64_t l = query_pos[i];

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, c);
    total_micros += cur_time_micros() - t0;
    total_ranks += ans;
  }

  std::cerr << "Total ranks: " << total_ranks << std::endl;
  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);
  write_log("Now running on an array of random positions for character G",
            LogLevel::MAJOR);

  total_ranks = 0;
  total_micros = 0;
  for (size_t i = 0; i < n; ++i) {
    if (query_pos[i] < 0) {
      throw std::runtime_error("Negative query length found");
    }
    if (query_pos[i] >= n_nodes) {
      throw std::runtime_error(
          "Query position out of bounds: " + std::to_string(query_pos[i]) +
          " >= " + std::to_string(n_nodes));
    }
    total_ranks %= 1000000007;

    char c = 'G';
    int64_t l = query_pos[i];

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, c);
    total_micros += cur_time_micros() - t0;
    total_ranks += ans;
  }

  std::cerr << "Total ranks: " << total_ranks << std::endl;
  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);
  write_log("Now running on an array of random positions for character T",
            LogLevel::MAJOR);
  total_micros = 0;
  for (size_t i = 0; i < n; ++i) {
    if (query_pos[i] < 0) {
      throw std::runtime_error("Negative query length found");
    }
    if (query_pos[i] >= n_nodes) {
      throw std::runtime_error(
          "Query position out of bounds: " + std::to_string(query_pos[i]) +
          " >= " + std::to_string(n_nodes));
    }
    total_ranks %= 1000000007;

    char c = 'T';
    int64_t l = query_pos[i];

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, c);
    total_micros += cur_time_micros() - t0;
    total_ranks += ans;
  }

  std::cerr << "Total ranks: " << total_ranks << std::endl;

  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);
  // Run
  std::sort(query_pos.begin(), query_pos.end());
  total_ranks = 0;
  total_micros = 0;
  for (size_t i = 0; i < n; ++i) {
    total_ranks %= 1000000007;

    char c = query_chars[i];
    int64_t l = query_pos[i];

    int64_t t0 = cur_time_micros();
    auto ans = s_rank_struct.rank(l, c);
    total_micros += cur_time_micros() - t0;
    total_ranks += ans;
    out_buffer.push_back(ans);
  }
  print_vector(out_buffer, writer);
  out_buffer.clear();
  std::cerr
      << "On an array of sorted positions (random characters), total ranks: "
      << total_ranks << std::endl;
  write_log("us/query: " + to_string((double)total_micros / number_of_queries) +
                " (excluding I/O etc)",
            LogLevel::MAJOR);
  return number_of_queries;
}

// Returns number of queries executed
template <typename sbwt_t>
int64_t run_queries(const string& infile_p, const string& infile_c,
                    const string& outfile, const sbwt_t& sbwt,
                    bool gzip_output) {
  typedef seq_io::Buffered_ofstream<seq_io::zstr::ofstream> out_gzip;
  typedef seq_io::Buffered_ofstream<std::ofstream> out_no_gzip;

  int64_t n_queries_run = 0;

  if (gzip_output) {
    n_queries_run +=
        run_file<sbwt_t, out_gzip>(infile_p, infile_c, outfile, sbwt);
  }
  if (!gzip_output) {
    n_queries_run +=
        run_file<sbwt_t, out_no_gzip>(infile_p, infile_c, outfile, sbwt);
  }

  return n_queries_run;
}

int s_rank_queries_main(int argc, char** argv) {
  int64_t micros_start = cur_time_micros();

  set_log_level(LogLevel::MINOR);

  cxxopts::Options options(
      argv[0],
      "Query just subset ranks at given positions for given characters.");

  options.add_options()("o,out-file", "Output filename.",
                        cxxopts::value<string>())(
      "i,index-file", "Index input file.", cxxopts::value<string>())(
      "p,pos-query-file", "The position query file.", cxxopts::value<string>())(
      "c,char-query-file", "The character query file.",
      cxxopts::value<string>())(
      "z,gzip-output",
      "Writes output in gzipped form. This can shrink the output files by an "
      "order of magnitude.",
      cxxopts::value<bool>()->default_value("false"))("h,help", "Print usage");

  int64_t old_argc = argc;  // Must store this because the parser modifies it
  auto opts = options.parse(argc, argv);

  if (old_argc == 1 || opts.count("help")) {
    std::cerr << options.help() << std::endl;
    exit(1);
  }

  string indexfile = opts["index-file"].as<string>();
  check_readable(indexfile);

  // Interpret input file
  string queryfile_c = opts["char-query-file"].as<string>();
  check_readable(queryfile_c);
  string queryfile_p = opts["pos-query-file"].as<string>();
  check_readable(queryfile_p);
  // Interpret output file
  string outfile = opts["out-file"].as<string>();
  bool gzip_output = opts["gzip-output"].as<bool>();
  check_writable(outfile);

  vector<string> variants = get_available_variants();

  throwing_ifstream in(indexfile, ios::binary);
  string variant = load_string(in.stream);  // read variant type
  if (std::find(variants.begin(), variants.end(), variant) == variants.end()) {
    cerr << "Error loading index from file: unrecognized variant specified in "
            "the file"
         << endl;
    return 1;
  }

  write_log("Loading the index variant " + variant, LogLevel::MAJOR);
  int64_t number_of_queries = 0;

  if (variant == "plain-matrix") {
    plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "rrr-matrix") {
    rrr_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "mef-matrix") {
    mef_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "plain-split") {
    plain_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "rrr-split") {
    rrr_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "mef-split") {
    mef_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "plain-concat") {
    plain_concat_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "mef-concat") {
    mef_concat_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "plain-subsetwt") {
    plain_sswt_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "rrr-subsetwt") {
    rrr_sswt_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "new-concat") {
    new_concat_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "correction-sets") {
    correction_sets_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "split-smaller-size") {
    split_smaller_size_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "concat-split-lengths") {
    concat_split_lengths_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "new-plain-concat") {
    new_plain_concat_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "new-split") {
    new_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "blocked-correction-sets") {
    blocked_correction_sets_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "fixed-block-correction-sets1") {
    fixed_block_correction_sets1_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "fixed-block-correction-sets2") {
    fixed_block_correction_sets2_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "fixed-block-correction-sets3") {
    fixed_block_correction_sets3_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "blocked8-split") {
    blocked8_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "blocked9-split") {
    blocked9_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  if (variant == "pred8-split") {
    pred8_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    number_of_queries +=
        run_queries(queryfile_p, queryfile_c, outfile, sbwt, gzip_output);
  }
  int64_t total_micros = cur_time_micros() - micros_start;
  write_log("us/query end-to-end: " +
                to_string((double)total_micros / number_of_queries),
            LogLevel::MAJOR);

  return 0;
}
