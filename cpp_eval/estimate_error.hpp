#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <stdexcept>

class BinnedSignedComplementaryGap{
private:
    std::vector<int> DB_values;
    std::discrete_distribution<int> distribution;
public:
    BinnedSignedComplementaryGap(const std::vector<std::pair<int,double>>& DB_freq_list) {
        if (DB_freq_list.empty()) {
            throw std::invalid_argument("empty sampling target");
        }
        std::vector<double> weights;
        for (const auto& DB_freq : DB_freq_list) {
            DB_values.push_back(DB_freq.first);
            weights.push_back(DB_freq.second);
        }
        distribution = std::discrete_distribution<int>(weights.begin(), weights.end());
    }

    int sampling_gap_in_DB(std::mt19937_64 &gen) {
        int DB_index = distribution(gen);
        int DB = DB_values[DB_index];
        return DB;
    }
    double sampling_error_rate(std::mt19937_64 &gen) {
        int DB = sampling_gap_in_DB(gen);
        double coefficient = 0.9; // confidence calibration
        double relative_error_rate = pow(10., -DB*coefficient/10.0);
        double error_rate = relative_error_rate / (1.0 + relative_error_rate);
        //std::cout << "DB: " << DB << ", relative_error_rate: " << relative_error_rate << ", error_rate: " << error_rate << std::endl;
        return error_rate;
    }
};

std::vector<double> create_log_probability_list(const std::vector<double>& error_rate_list) {
    int qubit_count = error_rate_list.size();
    std::vector<double> log_probability_list(1 << qubit_count);
    double sum = 0;
    for (int i = 0; i < (1 << qubit_count); ++i) {
        double log_probability = 0.0;
        for (int j = 0; j < qubit_count; ++j) {
            log_probability += log(1 - error_rate_list[j]);
            if ((i >> j) & 1) {
                log_probability += log(error_rate_list[j] / (1 - error_rate_list[j]));
            }
        }
        log_probability_list[i] = log_probability;
        sum += exp(log_probability);
    }
    return log_probability_list;
}

int calculate_syndrome(int error_pattern, int qubit_count) {
    int syndrome_bit = 0;
    for (int j = 0; j < qubit_count; ++j) {
        if ((error_pattern >> j) & 1) {
            syndrome_bit ^= 1;
        }
    }
    return syndrome_bit;
}

double optimal_decoder(std::vector<double> &error_rate_list) {
    int qubit_count = error_rate_list.size();

    // create a log probability list
    std::vector<double> log_probability_list = create_log_probability_list(error_rate_list);

    double logical_error_rate = 0.0;
    // check failure probability for each error pattern
    for (int error_pattern = 0; error_pattern < (1<<qubit_count); ++error_pattern) {

        // calculate syndrome of icerberg code
        int syndrome_bit = calculate_syndrome(error_pattern, qubit_count);

        // calculate the most probable class of error estimation
        double maximum_probability = 0.0;
        int most_probable_class = 0;
        for (int representative = 0; representative < (1<<(qubit_count-1)); ++representative) {
            // skip representatives that do not match the syndrome
            if (syndrome_bit != calculate_syndrome(representative, qubit_count)) {
                continue;
            }

            // update if the class with representative is the most probable one
            int opposite = ((1<<qubit_count)-1) ^ representative;
            double class_probability = exp(log_probability_list[representative]) + exp(log_probability_list[opposite]);
            if (class_probability > maximum_probability) {
                maximum_probability = class_probability;
                most_probable_class = representative;
            }
        }

        // accumulate if the most probable class does not match the error pattern or its complement
        double success = (most_probable_class == error_pattern) || (most_probable_class == (((1<<qubit_count)-1) ^ error_pattern));
        if (!success) {
            double error_probability = exp(log_probability_list[error_pattern]);
            logical_error_rate += error_probability;
        }
    }
    return logical_error_rate;
}

double minimum_weight_decoder(std::vector<double> &error_rate_list) {

    int qubit_count = error_rate_list.size();

    // create a log probability list
    std::vector<double> log_probability_list = create_log_probability_list(error_rate_list);

    double logical_error_rate = 0.0;
    // check failure probability for each error pattern
    for (int error_pattern = 0; error_pattern < (1<<qubit_count); ++error_pattern) {

        // calculate syndrome of icerberg code
        int syndrome_bit = calculate_syndrome(error_pattern, qubit_count);

        // calculate the most probable class of error estimation
        double maximum_probability = 0.0;
        int most_probable_error = 0;
        for (int error_estimation = 0; error_estimation < (1<<qubit_count); ++error_estimation) {
            // skip representatives that do not match the syndrome
            if (syndrome_bit != calculate_syndrome(error_estimation, qubit_count)) {
                continue;
            }

            // update if the error is the most probable one
            double error_probability = exp(log_probability_list[error_estimation]);
            if (error_probability > maximum_probability) {
                maximum_probability = error_probability;
                most_probable_error = error_estimation;
            }
        }

        // accumulate if the most probable class does not match the error pattern or its complement
        double success = (most_probable_error == error_pattern) || (most_probable_error == (((1<<qubit_count)-1) ^ error_pattern));
        if (!success) {
            double error_probability = exp(log_probability_list[error_pattern]);
            logical_error_rate += error_probability;
        }
    }
    return logical_error_rate;
}


std::pair<double, double> evaluate_logical_error_rate(int sample_count, int patch_count, int seed, BinnedSignedComplementaryGap &bscg, std::string policy, int num_rounds){
    if(patch_count % 2 == 1) {
        throw std::invalid_argument("patch_count must be even");
    }
    if(patch_count == 0) {
        throw std::invalid_argument("patch_count is zero");
    }

    std::mt19937_64 gen(seed);
    double sum_logical_error_rate = 0.0;
    std::vector<double> temp_logical_error_rate_list;
    for (int sample_index = 0; sample_index < sample_count; ++sample_index) {
        // sampling error rate list
        std::vector<double> error_rate_list(patch_count);
        for (int j = 0; j < patch_count; ++j) {
            error_rate_list[j] = bscg.sampling_error_rate(gen);
        }

        /*
        for (int j = 0; j < patch_count; ++j) {
            std::cout << error_rate_list[j] << " ";
        }
        std::cout << std::endl;
        */

        // evaluate failure probability based on the specified policy
        double temp_logical_error_rate = 0.0;
        if (policy == "optimal") {
            temp_logical_error_rate = optimal_decoder(error_rate_list);
        } else if (policy == "minimum_weight") {
            temp_logical_error_rate = minimum_weight_decoder(error_rate_list);
        } else {
            throw std::invalid_argument("Invalid policy: " + policy);
        }
        double temp_logical_error_rate_per_round = 1-pow(1-temp_logical_error_rate, 1./(num_rounds));
        double temp_logical_error_rate_per_patch_round = 1-pow(1-temp_logical_error_rate_per_round, 1./(patch_count));
        temp_logical_error_rate_per_patch_round *= 2; // account pL = px+pz

        temp_logical_error_rate_list.push_back(temp_logical_error_rate_per_patch_round);
        sum_logical_error_rate += temp_logical_error_rate_per_patch_round;
        // double current_sum_logical_error_rate = sum_logical_error_rate / (sample_index + 1);
        // std::cout << sample_index+1 << " " << temp_logical_error_rate << " " << current_sum_logical_error_rate << " +- " << std_error << std::endl;
    }
    double logical_error_rate = sum_logical_error_rate / sample_count;
    double std_error = 0.0;
    if (sample_count > 1) {
        double mean = logical_error_rate;
        double sum_squared_diff = 0.0;
        for (double x : temp_logical_error_rate_list) {
            double diff = x - mean;
            sum_squared_diff += diff * diff;
        }
        std_error = std::sqrt(sum_squared_diff / (sample_count - 1)) / std::sqrt(sample_count);
    }
    return std::make_pair(logical_error_rate, std_error);
}