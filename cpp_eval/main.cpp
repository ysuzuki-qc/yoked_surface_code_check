#include "estimate_error.hpp"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;


int main(int argc, char* argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <policy> <num_patch> <num_round_coef> <distance> <num_sample>" << endl;
        return 1;
    }

    std::string policy = argv[1];
    int num_patch = std::stoi(argv[2]);
    int num_round_coef = std::stoi(argv[3]);
    int distance = std::stoi(argv[4]);
    int num_sample = std::stoi(argv[5]);
    std::random_device rd;
    int seed = rd();

    std::string input_filename = "./data/gap_stat_" + std::to_string(num_round_coef) + "_" + std::to_string(distance) + ".txt";
    ifstream inputFile(input_filename);
    if (!inputFile) {
        cerr << "Error opening file: " << input_filename << endl;
        return 1;
    }

    std::vector<std::pair<int,double>> DB_freq_list;
    int DB;
    double freq;
    while (inputFile >> DB >> freq) {
        DB_freq_list.push_back({DB, freq});
        //std::cout << "Read DB: " << DB << ", freq: " << freq << std::endl;
    }

    BinnedSignedComplementaryGap complementary_gap_distribution(DB_freq_list);
    double logical_error_rate = evaluate_logical_error_rate(num_sample, num_patch, seed, complementary_gap_distribution, policy);
    double logical_error_rate_per_round = 1-pow(1-logical_error_rate, 1./(num_round_coef*distance));
    std::cout << "1D Yoke: n=" << num_patch << " r=" << num_round_coef << "d d=" << distance << ": logical error rate per round = " << logical_error_rate_per_round << std::endl;
    return 0;
}
