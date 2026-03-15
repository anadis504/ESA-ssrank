#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

static constexpr char alphabet[4] = {'A', 'C', 'G', 'T'};

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s [n] [output file] \n", argv[0]);
    std::exit(1);
  }
  int n = std::atoi(argv[1]);
  char* output_file = argv[2];
  std::vector<int64_t> indices;
  std::vector<char> letters;
  indices.reserve(1000000);
  letters.reserve(1000000);
  uint rng = 7;
  std::srand(rng);
  for (size_t i = 0; i < 1000000; ++i) {
    int64_t index = (int64_t)(std::rand() % n);
    char letter = alphabet[std::rand() % 4];
    indices.push_back(index);
    letters.push_back(letter);
  }

  auto ind_file = std::string(output_file) + "_indices.bin";
  auto letter_file = std::string(output_file) + "_letters.bin";
  std::ofstream ofs(ind_file, std::ios::binary);
  ofs.write(reinterpret_cast<char*>(indices.data()),
            sizeof(decltype(indices)::value_type) * indices.size());
  ofs.close();
  std::ofstream ofs2(letter_file, std::ios::binary);
  ofs2.write(reinterpret_cast<char*>(letters.data()),
             sizeof(decltype(letters)::value_type) * letters.size());
  ofs2.close();
}
