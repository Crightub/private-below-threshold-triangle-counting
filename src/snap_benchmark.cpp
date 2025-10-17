#include "counting.h"
#include "benchmark_utils.h"
#include "graph_generation.h"
#include "graph_statistics.h"
#include "logger.h"
#include "triangle.h"
#include "utils.h"


Graph load_snap_graph(const SnapBenchmarkConfig &cfg) {
    if (cfg.graph_name == "epinions") {
        return load_snap_epinions_graph("../snap/epinions.csv");
    } else if (cfg.graph_name == "wikipedia") {
        return load_snap_wikipedia_graph("../snap/wiki-rfa.txt");
    } else if (cfg.graph_name == "slashdot") {
        return load_snap_slashdot_graph("../snap/slashdot.txt");
    }

    return Graph(0);
}

void reset_loads(Graph &g) {
    auto load_map = get(&edge_info::load, g);

    for (auto e: boost::make_iterator_range(edges(g))) {
        put(load_map, e, 0);
    }
}

Benchmark run_snap_param_Benchmark(SnapBenchmarkConfig cfg) {
    Graph g = load_snap_graph(cfg);

    std::cout << "Loaded SNAP " << cfg.graph_name << " graph with " << boost::num_vertices(g) << " vertices and " <<
            boost::num_edges(g) <<
            " edges." << std::endl;

    // Graph structure is not altered so we find all triangles once in the beginning for a faster processing over multiple iterations
    std::list<Triangle> triangles = find_triangles(g);

    // Setup load balancing comparison -> only used if cfg.compare_load_balancing is activated
    std::map<double, Stat> results_load_balance;
    std::map<double, Stat> results_no_load_balance;

    auto load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, true, false};
    auto no_load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, false, true};

    for (int i = 0; i < cfg.iterations_count; ++i) {
        double value = cfg.start_value + i * cfg.step_size;
        set_param(cfg, cfg.param, value);

        std::cout << "Start iterations for " << to_string(cfg.param) << ": " << value << std::endl;
        Stat avg_load_balance;
        Stat avg_no_load_balance;

        for (int j = 0; j < cfg.instances_per_size; ++j) {
            Graph g_load_balance = add_discrete_laplace_noise(g, cfg.weight_eps);
            Graph g_no_load_balance = g_load_balance;
            reset_loads(g_load_balance);
            reset_loads(g_no_load_balance);

            std::cout << "Compute load balanced result..." << std::endl;
            auto triangles_load_balance = triangles;
            PrivateCountingResult res_load_balance = private_counting(g_load_balance, load_balancing_base_cfg,
                                                                      &triangles_load_balance);
            avg_load_balance.update(res_load_balance);

            if (cfg.compare_load_balancing) {
                std::cout << "Compute lowest index priority result..." << std::endl;
                auto triangles_no_load_balance = triangles;
                PrivateCountingResult res_no_load_balance = private_counting(
                    g_no_load_balance, no_load_balancing_base_cfg, &triangles_no_load_balance);
                avg_no_load_balance.update(res_no_load_balance);
            }
        }

        wrap_up_stats(results_load_balance, value, avg_load_balance);
        wrap_up_stats(results_no_load_balance, value, avg_no_load_balance);
    }

    return Benchmark{cfg, results_load_balance, results_no_load_balance};
}

int main(int argc, char *argv[]) {
    SnapBenchmarkConfig cfg = parse_snap_benchmark_config(argc, argv);

    std::string base_file_name = generate_filename(cfg);

    std::string log_filename = "../logs/" + base_file_name + ".log";
    std::string benchmark_filename = "../benchmark/" + base_file_name + ".json";

    init_logging(&log_filename);
    Benchmark bench = run_snap_param_Benchmark(cfg);

    save_benchmark_to_json(bench, benchmark_filename);

    return 0;
}
