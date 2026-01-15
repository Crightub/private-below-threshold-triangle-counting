#include "counting.h"
#include "benchmark_utils.h"
#include "graph_generation.h"
#include "graph_statistics.h"
#include "triangle.h"
#include "utils.h"


Graph load_graph(const BenchmarkConfig &cfg) {
    if (cfg.graph_name == "telecom-278") {
        return load_telecom_278_graph();
    }
    if (cfg.graph_name == "gmwcs") {
        return load_gmwcs_graph();
    }
    return {0};
}

void reset_loads(Graph &g) {
    auto load_map = get(&edge_info::load, g);

    for (auto e: boost::make_iterator_range(edges(g))) {
        put(load_map, e, 0);
    }
}


Benchmark run_instance_benchmark(BenchmarkConfig cfg) {
    Graph base_g = load_graph(cfg);

    std::cout << "Loaded SNAP " << cfg.graph_name << " graph with " << boost::num_vertices(base_g) << " vertices and "
            << boost::num_edges(base_g) << " edges." << std::endl;

    // Setup load balancing comparison -> only used if cfg.compare_load_balancing is activated
    std::map<double, Stat> results_load_balance;
    std::map<double, Stat> results_no_load_balance;

    auto load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, 4, true, false};
    auto no_load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, 4, false, true};

    for (int i = 0; i < cfg.iterations_count; ++i) {
        int value = cfg.start_value + i * cfg.step_size;
        set_param(load_balancing_base_cfg, cfg.param, value);
        set_param(no_load_balancing_base_cfg, cfg.param, value);

        Graph base_subgraph_g = extract_random_subgraph(base_g, value);

        auto g_load_balance = base_subgraph_g;
        auto g_no_load_balance = base_subgraph_g;

        std::cout << "Generated Subgraph with " << boost::num_vertices(base_subgraph_g) << " vertices and " <<
                boost::num_edges(base_subgraph_g) << " edges." << std::endl;

        std::vector<Triangle> triangles = find_triangles(base_subgraph_g);

        auto triangles_load_balance = triangles;
        auto triangles_no_load_balance = triangles;

        // Assign triangles based on load balancing
        for (Triangle &t: triangles_load_balance) {
            t.assign_triangle(g_load_balance, true);
        }

        if (cfg.compare_load_balancing) {
            for (Triangle &t: triangles_no_load_balance) {
                t.assign_triangle(g_no_load_balance, false);
            }
        }

        std::cout << "Start iterations for " << to_string(cfg.param) << ": " << value << std::endl;
        auto avg_load_balance = Stat(cfg.instances_per_iter);
        auto avg_no_load_balance = Stat(cfg.instances_per_iter);

        for (int j = 0; j < cfg.instances_per_iter; ++j) {
            add_discrete_laplace_noise(g_load_balance, cfg.weight_eps, cfg.count_eps);
            std::cout << "Iteration: " << j << std::endl;
            std::cout << "Compute load balanced result..." << std::endl;
            PrivateCountingResult res_load_balance = private_counting(g_load_balance, load_balancing_base_cfg,
                                                                      &triangles_load_balance);
            avg_load_balance.update(res_load_balance);

            if (cfg.compare_load_balancing) {
                std::cout << "Compute lowest index priority result..." << std::endl;
                add_discrete_laplace_noise(g_no_load_balance, cfg.weight_eps, cfg.count_eps);

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

Benchmark run_param_benchmark(BenchmarkConfig cfg) {
    Graph g_load_balance = load_graph(cfg);
    Graph g_no_load_balance = g_load_balance;

    std::cout << "Loaded " << cfg.graph_name << " graph with " << boost::num_vertices(g_load_balance) <<
            " vertices and " <<
            boost::num_edges(g_load_balance) <<
            " edges." << std::endl;

    // Graph structure is not altered so we find all triangles once in the beginning for a faster processing over multiple iterations
    std::vector<Triangle> triangles_load_balance = find_triangles(g_load_balance);
    std::vector<Triangle> triangles_no_load_balance = triangles_load_balance;

    std::cout << "Assign triangles to nodes responsible for counting." << std::endl;

    // Assign triangles based on load balancing
    for (Triangle &t: triangles_load_balance) {
        t.assign_triangle(g_load_balance, true);
    }

    if (cfg.compare_load_balancing) {
        for (Triangle &t: triangles_no_load_balance) {
            t.assign_triangle(g_no_load_balance, false);
        }
    }

    // Setup load balancing comparison -> only used if cfg.compare_load_balancing is activated
    std::map<double, Stat> results_load_balance;
    std::map<double, Stat> results_no_load_balance;

    auto load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, 4, true};
    auto no_load_balancing_base_cfg = PrivateCountingConfig{cfg.weight_eps, cfg.count_eps, cfg.lambda, 4, false};

    for (int i = 0; i < cfg.iterations_count; ++i) {
        double value = cfg.start_value + i * cfg.step_size;
        set_param(load_balancing_base_cfg, cfg.param, value);
        set_param(no_load_balancing_base_cfg, cfg.param, value);

        auto avg_load_balance = Stat(cfg.instances_per_iter);
        auto avg_no_load_balance = Stat(cfg.instances_per_iter);

        for (int j = 0; j < cfg.instances_per_iter; ++j) {
            std::cout << to_string(cfg.param) << ": " << value << ", Iteration: " << j << std::endl;
            add_discrete_laplace_noise(g_load_balance, cfg.weight_eps, cfg.count_eps);

            std::cout << "Compute load balanced result..." << std::endl;
            PrivateCountingResult res_load_balance = private_counting(g_load_balance,
                                                                      load_balancing_base_cfg,
                                                                      &triangles_load_balance);
            avg_load_balance.update(res_load_balance);

            if (cfg.compare_load_balancing) {
                std::cout << "Compute lowest index priority result..." << std::endl;
                add_discrete_laplace_noise(g_no_load_balance, cfg.weight_eps, cfg.count_eps);
                PrivateCountingResult res_no_load_balance = private_counting(
                    g_no_load_balance,
                    no_load_balancing_base_cfg,
                    &triangles_no_load_balance);
                avg_no_load_balance.update(res_no_load_balance);
            }
        }

        wrap_up_stats(results_load_balance, value, avg_load_balance);
        wrap_up_stats(results_no_load_balance, value, avg_no_load_balance);
    }

    return Benchmark{cfg, results_load_balance, results_no_load_balance};
}

int main(int argc, char *argv[]) {
    BenchmarkConfig cfg = parse_benchmark_config(argc, argv);
    cfg.compare_load_balancing = true;
    std::string base_file_name = generate_filename(cfg);

    std::string benchmark_filename = "../benchmark/" + base_file_name + ".json";

    Benchmark bench;
    if (cfg.param == Param::InstanceSize) {
        std::cout << "Run instance benchmark" << std::endl;
        bench = run_instance_benchmark(cfg);
    } else {
        bench = run_param_benchmark(cfg);
    }

    save_benchmark_to_json(bench, benchmark_filename);

    return 0;
}
