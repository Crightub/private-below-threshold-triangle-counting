#pragma once

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

enum class Param {
    InstanceSize,
    WeightEps,
    CountEps,
    Lambda,
    Gamma
};

inline std::string to_string(Param p) {
    switch (p) {
        case Param::InstanceSize:
            return "instance_size";
        case Param::WeightEps:
            return "weight_eps";
        case Param::CountEps:
            return "count_eps";
        case Param::Lambda:
            return "lambda";
        case Param::Gamma:
            return "gamma";
        default:
            return "unknown";
    }
}

const std::unordered_map<std::string, Param> stringToParam = {
    {"instance_size", Param::InstanceSize},
    {"weight_eps", Param::WeightEps},
    {"count_eps", Param::CountEps},
    {"lambda", Param::Lambda},
    {"gamma", Param::Gamma}
};

struct Stat {
    long opt = 0;

    double naive = 0;
    double naive_l2 = 0;

    double global_unbiased = 0;
    double global_biased = 0;
    double smooth_unbiased = 0;
    double smooth_biased = 0;

    double global_unbiased_l2 = 0;
    double global_biased_l2 = 0;
    double smooth_unbiased_l2 = 0;
    double smooth_biased_l2 = 0;

    int n = 0;

    std::vector<PrivateCountingResult> results;

    Stat() = default;

    explicit Stat(int _n) : n(_n) {
    }

    void update(const PrivateCountingResult &res) {
        opt = res.opt;
        naive += 1.0 / n * res.naive;
        global_unbiased += 1.0 / n * res.global_unbiased;
        global_biased += 1.0 / n * res.global_biased;
        smooth_unbiased += 1.0 / n * res.smooth_unbiased;
        smooth_biased += 1.0 / n * res.smooth_biased;

        naive_l2 += 1.0 / n * (res.naive - res.opt) * (res.naive - res.opt);
        global_unbiased_l2 += 1.0 / n * (res.global_unbiased - res.opt) * (res.global_unbiased - res.opt);
        global_biased_l2 += 1.0 / n * (res.global_biased - res.opt) * (res.global_biased - res.opt);
        smooth_unbiased_l2 += 1.0 / n * (res.smooth_unbiased - res.opt) * (res.smooth_unbiased - res.opt);
        smooth_biased_l2 += 1.0 / n * (res.smooth_biased - res.opt) * (res.smooth_biased - res.opt);

        results.push_back(res);
    }
};


struct BenchmarkConfig {
    std::string graph_name;
    int instances_per_iter;
    int iterations_count;
    double start_value;
    double step_size;
    Param param;
    bool compare_load_balancing;

    int instance_size;
    double weight_eps;
    double count_eps;
    int lambda;
};

inline void print_benchmark_config(const BenchmarkConfig &cfg) {
    std::cout << "Benchmark Config\n------------------\n"
            << "Graph: " << cfg.graph_name << std::endl
            << "Flexible Parameter: " << to_string(cfg.param) << std::endl
            << "Instances per Iteration: " << cfg.instances_per_iter << std::endl
            << "Iterations: " << cfg.iterations_count << std::endl
            << "Step Size: " << cfg.step_size << std::endl
            << "Weight eps: " << cfg.weight_eps << std::endl
            << "Count eps: " << cfg.count_eps << std::endl
            << "Lambda: " << cfg.lambda << std::endl
            << "------------------" << std::endl;
}

struct Benchmark {
    BenchmarkConfig cfg;
    std::map<double, Stat> results_load_balance;
    std::map<double, Stat> results_no_load_balance;
};

inline void to_json(json &j, const BenchmarkConfig &cfg) {
    j = {
        {"graph_name", cfg.graph_name},
        {"parameter", to_string(cfg.param)},
        {"instances_per_iter", cfg.instances_per_iter},
        {"iterations_count", cfg.iterations_count},
        {"step_size", cfg.step_size},
        {"weight_eps", cfg.weight_eps},
        {"count_eps", cfg.count_eps},
        {"lambda", cfg.lambda}
    };
}

inline void to_json(json &j, const Stat &stat) {
    j = {
        {"opt", stat.opt},
        {"global_unbiased", stat.global_unbiased},
        {"global_biased", stat.global_biased},
        {"smooth_unbiased", stat.smooth_unbiased},
        {"smooth_biased", stat.smooth_biased},
        {"global_unbiased_l2", stat.global_unbiased_l2},
        {"global_biased_l2", stat.global_biased_l2},
        {"smooth_unbiased_l2", stat.smooth_unbiased_l2},
        {"smooth_biased_l2", stat.smooth_biased_l2},
        {"naive", stat.naive},
        {"naive_l2", stat.naive_l2},
        {"results", stat.results}
    };
}

inline BenchmarkConfig parse_benchmark_config(int argc, char *argv[]) {
    if (argc != 10 && argc != 6) {
        std::cerr << "Single Iteration Benchmark: " << argv[0] <<
                " <graph_name> <instances> <weight_eps> <count_eps> <lambda>" << std::endl
                << "Multi Iteration Benchmark: " << argv[0]
                << " <graph_name> <instances_per_size> <iterations_count> <start_value> <step_size> <param> <weight_eps> <count_eps> <lambda>"
                << std::endl;
        exit(1);
    }

    if (argc == 10) {
        BenchmarkConfig settings{};

        settings.graph_name = argv[1];
        settings.instances_per_iter = std::atoi(argv[2]);
        settings.iterations_count = std::atoi(argv[3]);
        settings.start_value = std::atof(argv[4]);
        settings.step_size = std::atof(argv[5]);
        settings.param = stringToParam.at(argv[6]);
        settings.weight_eps = std::atof(argv[7]);
        settings.count_eps = std::atof(argv[8]);
        settings.lambda = std::atoi(argv[9]);
        return settings;
    }


    BenchmarkConfig settings{};

    settings.graph_name = argv[1];
    settings.instances_per_iter = std::atoi(argv[2]);
    settings.weight_eps = std::atof(argv[3]);
    settings.count_eps = std::atof(argv[4]);
    settings.lambda = std::atoi(argv[5]);

    // Set dummy setting
    settings.iterations_count = 1;
    settings.start_value = settings.weight_eps;
    settings.step_size = 0;
    settings.param = Param::WeightEps;

    return settings;
}

inline void set_param(PrivateCountingConfig &cfg, const Param param, double value) {
    switch (param) {
        case Param::WeightEps:
            cfg.weight_eps = value;
            break;
        case Param::CountEps:
            cfg.count_eps = value;
            break;
        case Param::Lambda:
            cfg.lambda = value;
            break;
        case Param::Gamma:
            cfg.gamma = value;
        default:
            break;
    }
}


inline void wrap_up_stats(std::map<double, Stat> &results, double value, Stat avg) {
    results[value] = avg;
}


inline void save_benchmark_to_json(const Benchmark &bench, const std::string &filename) {
    std::cout << "Save benchmark results at " << filename << std::endl;
    json j;
    j["config"] = bench.cfg;

    for (const auto &[key, result]: bench.results_load_balance) {
        j["results_map"]["load_balance"][std::to_string(key)] = result;
    }

    for (const auto &[key, result]: bench.results_no_load_balance) {
        j["results_map"]["no_load_balance"][std::to_string(key)] = result;
    }

    std::ofstream out(filename);
    out << j.dump(4);
}


inline std::string generate_filename(const BenchmarkConfig &cfg) {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);

    std::ostringstream timestamp;
    timestamp << std::put_time(&tm, "%Y%m%d_%H%M");

    std::ostringstream filename;
    filename << "bench_"
            << cfg.graph_name << "_"
            << to_string(cfg.param) << "_"
            << cfg.start_value << "_"
            << cfg.iterations_count << "_"
            << cfg.step_size << "_"
            << cfg.weight_eps << "_"
            << cfg.count_eps << "_"
            << cfg.lambda << "_"
            << timestamp.str();

    return filename.str();
}

inline Graph extract_random_subgraph(const Graph &graph, std::size_t size) {
    const std::size_t n = boost::num_vertices(graph);
    if (size > n) {
        throw std::invalid_argument("Requested subgraph larger than original graph");
    }

    std::vector<Node> vertices;
    vertices.reserve(n);
    for (Node v = 0; v < n; ++v)
        vertices.push_back(v);

    std::shuffle(vertices.begin(), vertices.end(), rng);
    vertices.resize(size);

    std::unordered_map<Node, std::size_t> index_map;
    for (std::size_t i = 0; i < size; ++i)
        index_map[vertices[i]] = i;

    Graph subgraph(size);

    // Copy edge weights
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = i + 1; j < size; ++j) {
            auto e_pair = boost::edge(vertices[i], vertices[j], graph);
            if (e_pair.second) {
                const edge_info &info = graph[e_pair.first];
                boost::add_edge(index_map[vertices[i]], index_map[vertices[j]], info, subgraph);
            }
        }
    }

    return subgraph;
}
