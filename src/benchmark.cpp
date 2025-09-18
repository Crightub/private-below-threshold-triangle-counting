#include <iostream>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include "utils.h"
#include "distribution.h"
#include "graph_generation.h"
#include "qp_counting.h"
#include "basic_counting.h"
#include <map>
#include <unordered_map>
#include <nlohmann/json.hpp>

#include "logger.h"
#include "../test/graph_presets.h"
using json = nlohmann::json;

enum class Param {
    InstanceSize,
    EdgeProbability,
    WeightMu,
    WeightStd,
    WeightEps,
    CountEps
};

std::string to_string(Param p) {
    switch (p) {
        case Param::InstanceSize:
            return "instance_size";
        case Param::EdgeProbability:
            return "edge_probability";
        case Param::WeightMu:
            return "weight_mu";
        case Param::WeightStd:
            return "weight_std";
        case Param::WeightEps:
            return "weight_eps";
        case Param::CountEps:
            return "count_eps";
        default:
            return "unknown";
    }
}

const std::unordered_map<std::string, Param> stringToParam = {
    {"instance_size", Param::InstanceSize},
    {"edge_probability", Param::EdgeProbability},
    {"weight_mu", Param::WeightMu},
    {"weight_std", Param::WeightStd},
    {"weight_eps", Param::WeightEps},
    {"count_eps", Param::CountEps}
};

struct Stat {
    double count_avg = 0.0;
    double weight_noise = 0.0;
    double count_noise = 0.0;

    int n = 0;

    void update(const PrivateCountingResult &res) {
        count_avg += res.negative_triangle_count;
        weight_noise += res.weight_noise;
        count_noise += res.count_noise;
        n++;
    }

    void normalize() {
        count_avg /= n;
        weight_noise /= n;
        count_noise /= n;
    }
};


struct BenchmarkResult {
    double opt_negative_triangle_count;
    Stat unbiased_global;
    Stat unbiased_smooth;
    Stat biased_global;
    Stat biased_smooth;
};

struct BenchmarkConfig {
    int instances_per_size;
    int iterations_count;
    double start_value;
    double step_size;
    Param param;

    int instance_size;
    double edge_probability;
    double weight_mu;
    double weight_std;
    double weight_eps;
    double count_eps;
};

void print_benchmark_config(const BenchmarkConfig &cfg) {
    std::cout << "Benchmark Config\n------------------\n"
            << "Flexible Parameter: " << to_string(cfg.param) << std::endl
            << "Instances per Size: " << cfg.instances_per_size << std::endl
            << "Iterations: " << cfg.iterations_count << std::endl
            << "Step Size: " << cfg.step_size << std::endl
            << "Edge Probability: " << cfg.edge_probability << std::endl
            << "Weight Mu: " << cfg.weight_mu << std::endl
            << "Weight Std: " << cfg.weight_std << std::endl
            << "Weight eps: " << cfg.weight_eps << std::endl
            << "Count eps: " << cfg.count_eps << std::endl
            << "------------------" << std::endl;
}

struct Benchmark {
    BenchmarkConfig cfg;
    std::map<double, BenchmarkResult> results_map;
};

void to_json(json &j, const BenchmarkConfig &cfg) {
    j = {
        {"parameter", to_string(cfg.param)},
        {"instances_per_size", cfg.instances_per_size},
        {"iterations_count", cfg.iterations_count},
        {"step_size", cfg.step_size},
        {"edge_probability", cfg.edge_probability},
        {"weight_mu", cfg.weight_mu},
        {"weight_std", cfg.weight_std},
        {"weight_eps", cfg.weight_eps},
        {"count_eps", cfg.count_eps}
    };
}

void to_json(json &j, const Stat &stat) {
    j = {
        {"count_avg", stat.count_avg},
        {"weight_noise", stat.weight_noise},
        {"count_noise", stat.count_noise}
    };
}

void to_json(json &j, const BenchmarkResult &res) {
    j = {
        {"opt_negative_triangle_count", res.opt_negative_triangle_count},
        {"unbiased_global", res.unbiased_global},
        {"unbiased_smooth", res.unbiased_smooth},
        {"biased_global", res.biased_global},
        {"biased_smooth", res.biased_smooth}
    };
}


BenchmarkConfig parse_benchmark_config(int argc, char *argv[]) {
    if (argc < 12) {
        std::cerr << "Usage: " << argv[0]
                << " <instances_per_size> <iterations_count> <start_value> <step_size> <param> <instance_size> "
                << "<edge_probability> <weight_mu> <weight_std> "
                << "<weight_eps> <count_eps>" << std::endl;
        std::exit(1);
    }

    BenchmarkConfig settings{};

    settings.instances_per_size = std::atoi(argv[1]);
    settings.iterations_count = std::atoi(argv[2]);
    settings.start_value = std::atof(argv[3]);
    settings.step_size = std::atof(argv[4]);
    settings.param = stringToParam.at(argv[5]);
    settings.instance_size = std::atoi(argv[6]);
    settings.edge_probability = std::atof(argv[7]);
    settings.weight_mu = std::atof(argv[8]);
    settings.weight_std = std::atof(argv[9]);
    settings.weight_eps = std::atof(argv[10]);
    settings.count_eps = std::atof(argv[11]);

    return settings;
}

void set_param(BenchmarkConfig &cfg, const Param param, double value) {
    switch (param) {
        case Param::InstanceSize:
            cfg.instance_size = value;
            break;
        case Param::EdgeProbability:
            cfg.edge_probability = value;
            break;
        case Param::WeightMu:
            cfg.weight_mu = value;
            break;
        case Param::WeightStd:
            cfg.weight_std = value;
            break;
        case Param::WeightEps:
            cfg.weight_eps = value;
            break;
        case Param::CountEps:
            cfg.count_eps = value;
            break;
        default:
            throw std::invalid_argument("Unknown parameter name.");
    }
}


void wrap_up_stats(std::map<double, BenchmarkResult> &results, double value, Stat unbiased_global, Stat unbiased_smooth,
                   Stat biased_global, Stat biased_smooth, double opt_avg) {
    unbiased_global.normalize();
    unbiased_smooth.normalize();
    biased_global.normalize();
    biased_smooth.normalize();

    results[value] = BenchmarkResult{
        opt_avg,
        unbiased_global,
        unbiased_smooth,
        biased_global,
        biased_smooth
    };
}

Graph generate_preset_graph(int instance_size, double weight_mu, double weight_std, double weight_eps) {
    Graph g = star_outline_graph(instance_size);
    g = assign_gaussian_weights(g, weight_mu, weight_std);
    Graph noisyGraph = add_discrete_laplace_noise(g, weight_eps);

    return noisyGraph;
}

Benchmark run_param_benchmark(BenchmarkConfig cfg) {
    std::map<double, BenchmarkResult> results;

    print_benchmark_config(cfg);

    const PrivateCountingConfig unbiased_global_cfg{true, cfg.weight_eps, cfg.count_eps, false};
    const PrivateCountingConfig unbiased_smooth_cfg{true, cfg.weight_eps, cfg.count_eps, true};
    const PrivateCountingConfig biased_global_cfg{false, cfg.weight_eps, cfg.count_eps, false};
    const PrivateCountingConfig biased_smooth_cfg{false, cfg.weight_eps, cfg.count_eps, true};

    for (int i = 0; i < cfg.iterations_count; ++i) {
        double value = cfg.start_value + i * cfg.step_size;
        set_param(cfg, cfg.param, value);

        Stat unbiased_global_stat;
        Stat unbiased_smooth_stat;
        Stat biased_global_stat;
        Stat biased_smooth_stat;

        double opt_avg = 0.0;

        std::cout << "Start instance loop..." << std::endl;
        // #pragma omp parallel for
        for (int j = 0; j < cfg.instances_per_size; ++j) {
            Graph g = generate_graph(cfg.instance_size, cfg.edge_probability, cfg.weight_mu, cfg.weight_std,
                                     cfg.weight_eps);
            // Graph g = generate_preset_graph(cfg.instance_size, cfg.weight_mu, cfg.weight_std, cfg.weight_eps);
            // std::cout << "Generated graph for param = " << value << std::endl;

            int opt_negative_triangles = count_negative_triangles(g);
            std::cout << "opt negative triangles: " << opt_negative_triangles << std::endl;
            PrivateCountingResult unbiased_global_result = randomized_private_counting(g, unbiased_global_cfg);
            std::cout << "unbiased_global: " << unbiased_global_result.negative_triangle_count << std::endl;
            // PrivateCountingResult unbiased_smooth_result = randomized_private_counting(g, unbiased_smooth_cfg);
            // std::cout << "unbiased_smooth: " << unbiased_smooth_result.negative_triangle_count << std::endl;
            PrivateCountingResult biased_global_result = randomized_private_counting(g, biased_global_cfg);
            std::cout << "biased_global: " << biased_global_result.negative_triangle_count << std::endl;
            // PrivateCountingResult biased_smooth_result = randomized_private_counting(g, biased_smooth_cfg);
            // std::cout << "biased_smooth: " << biased_smooth_result.negative_triangle_count << std::endl;

            opt_avg += opt_negative_triangles;
            unbiased_global_stat.update(unbiased_global_result);
            // unbiased_smooth_stat.update(unbiased_smooth_result);
            biased_global_stat.update(biased_global_result);
            // biased_smooth_stat.update(biased_smooth_result);
        }

        opt_avg /= cfg.instances_per_size;
        wrap_up_stats(results, value, unbiased_global_stat, unbiased_smooth_stat, biased_global_stat,
                      biased_smooth_stat, opt_avg);
    }

    return Benchmark{cfg, results};
}

void save_benchmark_to_json(const Benchmark &bench, const std::string &filename) {
    std::cout << "Save benchmark results at " << filename << std::endl;
    json j;
    j["config"] = bench.cfg;

    for (const auto &[key, result]: bench.results_map) {
        j["results_map"][std::to_string(key)] = result;
    }

    std::ofstream out(filename);
    out << j.dump(4);
}

std::string generate_filename(const BenchmarkConfig &cfg) {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);

    std::ostringstream timestamp;
    timestamp << std::put_time(&tm, "%Y%m%d_%H%M");

    std::ostringstream filename;
    filename << "bench_"
            << to_string(cfg.param) << "_"
            << cfg.iterations_count << "_"
            << cfg.step_size << "_"
            << cfg.instance_size << "_"
            << cfg.edge_probability << "_"
            << cfg.weight_mu << "_"
            << cfg.weight_eps << "_"
            << cfg.count_eps << "_"
            << timestamp.str();

    return filename.str();
}

int main(int argc, char *argv[]) {
    BenchmarkConfig cfg = parse_benchmark_config(argc, argv);
    std::string base_file_name = generate_filename(cfg);

    std::string log_filename = "../logs/" + base_file_name + ".log";
    std::string benchmark_filename = "../benchmark/" + base_file_name + ".json";

    init_logging(&log_filename);
    Benchmark bench = run_param_benchmark(cfg);

    save_benchmark_to_json(bench, benchmark_filename);

    return 0;
}
