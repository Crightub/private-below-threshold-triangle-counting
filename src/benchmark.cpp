#include <boost/graph/graphviz.hpp>
#include "utils.h"
#include "graph_generation.h"
#include "counting.h"
#include <map>
#include <nlohmann/json.hpp>

#include "benchmark_utils.h"
#include "graph_statistics.h"
#include "logger.h"
#include "../test/graph_presets.h"
using json = nlohmann::json;

//
// Graph generate_preset_graph(int instance_size, double weight_mu, double weight_std, double weight_eps) {
//     Graph g = complete_graph(instance_size);
//     g = assign_gaussian_weights(g, weight_mu, weight_std);
//     Graph noisyGraph = add_discrete_laplace_noise(g, weight_eps);
//
//     return noisyGraph;
// }
//
// Benchmark run_param_benchmark(BenchmarkConfig cfg) {
//     std::map<double, Stat> results_load_balance;
//     std::map<double, Stat> results_no_load_balance;
//
//     print_benchmark_config(cfg);
//
//     auto load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, true, false};
//     auto no_load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, false, true};
//
//     for (int i = 0; i < cfg.iterations_count; ++i) {
//         double value = cfg.start_value + i * cfg.step_size;
//         set_param(cfg, cfg.param, value);
//
//         Stat load_balance_avg;
//         Stat no_load_balance_avg;
//
//         for (int j = 0; j < cfg.instances_per_size; ++j) {
//             std::cout << "Generate graph for " << to_string(cfg.param) << " = " << value << std::endl;
//             Graph g_load_balance = generate_graph(cfg.instance_size, cfg.edge_probability, cfg.weight_mu,
//                                                   cfg.weight_std,
//                                                   cfg.weight_eps);
//             Graph g_no_load_balance = g_load_balance;
//
//             std::list<Triangle> triangles = find_triangles(g_load_balance);
//
//             Graph g = generate_preset_graph(cfg.instance_size, cfg.weight_mu, cfg.weight_std, cfg.weight_eps);
//             std::cout << "Compute load balanced result..." << std::endl;
//             PrivateCountingResult res_load_balance = private_counting(g_load_balance, load_balancing_base_cfg,
//                                                                       &triangles);
//             load_balance_avg.update(res_load_balance);
//
//             std::cout << "Compute not load balanced result..." << std::endl;
//             PrivateCountingResult res_no_load_balance = private_counting(
//                 g_no_load_balance, no_load_balancing_base_cfg, &triangles);
//             no_load_balance_avg.update(res_no_load_balance);
//         }
//
//         wrap_up_stats(results_load_balance, value, load_balance_avg);
//         wrap_up_stats(results_no_load_balance, value, no_load_balance_avg);
//     }
//
//     return Benchmark{cfg, results_load_balance, results_no_load_balance};
// }
//
//
// int main(int argc, char *argv[]) {
//     BenchmarkConfig cfg = parse_benchmark_config(argc, argv);
//     std::string base_file_name = generate_filename(cfg);
//
//     std::string log_filename = "../logs/" + base_file_name + ".log";
//     std::string benchmark_filename = "../benchmark/" + base_file_name + ".json";
//
//     init_logging(&log_filename);
//     Benchmark bench = run_param_benchmark(cfg);
//
//     save_benchmark_to_json(bench, benchmark_filename);
//
//     return 0;
// }
