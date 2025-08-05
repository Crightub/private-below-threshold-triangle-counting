#include <iostream>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include "utils.hpp"
#include "distribution.hpp"
#include "graph_generation.hpp"
#include "qp_counting.hpp"
#include "basic_counting.hpp"
#include <cmath>
#include <map>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

struct BenchmarkResult
{
    double opt_negative_triangle_count;
    double unbiased_private_triangle_count;
    double biased_private_triangle_count;
};

struct BenchmarkConfig
{
    int instances_per_size;
    int max_instance_size;
    int step_size;

    double edge_probability;
    double weight_mu;
    double weight_std;
    double weight_eps;
    double count_eps;
};

void print_benchmark_config(const BenchmarkConfig& cfg){
    std::cout << "Benchmark Config\n------------------\n" 
    << "Instances per Size: " << cfg.instances_per_size << std::endl
    << "Max Instance Size: " << cfg.max_instance_size << std::endl
    << "Step Size: " << cfg.step_size << std::endl
    << "Edge Probability: " << cfg.edge_probability << std::endl
    << "Weight Mu: " << cfg.weight_mu << std::endl
    << "Weight Std: " << cfg.weight_std << std::endl
    << "Weight eps: " <<  cfg.weight_eps << std::endl
    << "Count eps: " << cfg.count_eps << std::endl
    << "------------------" << std::endl;
}

struct Benchmark
{
    BenchmarkConfig cfg;
    std::map<int, BenchmarkResult> results_map;
};

void to_json(json& j, const BenchmarkConfig& cfg)
{
    j = {
        {"instances_per_size", cfg.instances_per_size},
        {"max_instance_size", cfg.max_instance_size},
        {"step_size", cfg.step_size},
        {"edge_probability", cfg.edge_probability},
        {"weight_mu", cfg.weight_mu},
        {"weight_std", cfg.weight_std},
        {"weight_eps", cfg.weight_eps},
        {"count_eps", cfg.count_eps}
    };
}

void to_json(json& j, const BenchmarkResult& res)
{
    j = {
        {"opt_negative_triangle_count", res.opt_negative_triangle_count},
        {"unbiased_private_triangle_count", res.unbiased_private_triangle_count},
        {"biased_private_triangle_count", res.biased_private_triangle_count}
    };
}


BenchmarkConfig parse_benchmark_config(int argc, char* argv[])
{
    if (argc < 9)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <instances_per_size> <max_instance_size> <step_size>"
                  << " <edge_probability> <weight_mu> <weight_std>"
                  << " <weight_eps> <count_eps>" << std::endl;
        std::exit(1);
    }

    BenchmarkConfig settings;

    settings.instances_per_size = std::atoi(argv[1]);
    settings.max_instance_size = std::atoi(argv[2]);
    settings.step_size = std::atoi(argv[3]);

    settings.edge_probability = std::atof(argv[4]);
    settings.weight_mu = std::atof(argv[5]);
    settings.weight_std = std::atof(argv[6]);
    settings.weight_eps = std::atof(argv[7]);
    settings.count_eps = std::atof(argv[8]);

    return settings;
}

Benchmark run_benchmark(BenchmarkConfig cfg){

    std::map<int, BenchmarkResult> result = {};

    print_benchmark_config(cfg);

    for (int i = cfg.step_size; i < cfg.max_instance_size; i += cfg.step_size)
    {
        double avg_opt_count = 0;
        double avg_biased_count = 0;
        double avg_unbiased_count = 0;

        for (int j = 0; j < cfg.instances_per_size; ++j)
        {
            Graph g = generateGraph(i, cfg.edge_probability, cfg.weight_mu, cfg.weight_std, cfg.weight_eps);

            int opt_negative_triangles = count_negative_triangles(g);
            PrivateCountingResult unbiased_result = randomized_private_counting(g, cfg.weight_eps, cfg.count_eps, true, true, false);
            PrivateCountingResult biased_result = randomized_private_counting(g, cfg.weight_eps, cfg.count_eps, true, true, true);

            avg_opt_count += opt_negative_triangles;
            avg_unbiased_count = unbiased_result.negative_triangle_count;
            avg_biased_count = biased_result.negative_triangle_count;
        }

        avg_opt_count /= cfg.instances_per_size;
        avg_unbiased_count /= cfg.instances_per_size;
        avg_biased_count /= cfg.instances_per_size;

        std::cout << "n = " << i << ", opt: " << avg_opt_count << ", unbiased: " << avg_unbiased_count << ", biased: " << avg_biased_count << std::endl;
        
        result[i] = BenchmarkResult{avg_opt_count, avg_unbiased_count, avg_biased_count};
    }

    Benchmark bench_out = Benchmark{
        cfg,
        result
    };

    return bench_out;
}

void save_benchmark_to_json(const Benchmark& bench, const std::string& filename)
{
    std::cout << "Save benchmark results at " << filename << std::endl;
    json j;
    j["config"] = bench.cfg;

    for (const auto& [key, result] : bench.results_map)
    {
        j["results_map"][std::to_string(key)] = result;
    }

    std::ofstream out(filename);
    out << j.dump(4);
}

std::string generate_filename(const BenchmarkConfig& cfg)
{
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);

    std::ostringstream timestamp;
    timestamp << std::put_time(&tm, "%Y%m%d_%H%M");

    std::ostringstream filename;
    filename << "bench_"
             << cfg.step_size << "_"
             << cfg.max_instance_size << "_"
             << cfg.edge_probability << "_"
             << cfg.weight_mu << "_"
             << cfg.weight_eps << "_"
             << cfg.count_eps << "_"
             << timestamp.str()
             << ".json";

    return filename.str();
}


int main(int argc, char *argv[])
{
    BenchmarkConfig cfg = parse_benchmark_config(argc, argv);

    Benchmark bench = run_benchmark(cfg);

    std::string filename = generate_filename(cfg);

    save_benchmark_to_json(bench, "../benchmark/" + filename);

    return 0;
}
