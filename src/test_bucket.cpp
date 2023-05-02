#include <random>
#include <chrono>
#include <iostream>
#include "ParticleLeniaBucket.h"
#include "ParticleLeniaBucketThreaded.h"

int main() {
    int num_particles = 2;
    bool human_readable = false;
    long max_time = 120 * 1000;

    while (true) {
        if (human_readable) {
            std::cout << "=== Evaluating n=" << num_particles << " case ===" << std::endl;
        }
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_real_distribution<> distribution(0, 1.1 * std::sqrt(num_particles));

        std::vector<Point> particles;
        for (int i = 0; i < num_particles; ++i)
            particles.push_back(Point{
                    distribution(rng), distribution(rng)
            });

        Parameters parameters;
        //ParticleLeniaBucket particle_lenia(parameters, particles, 1e-4, 0.5, human_readable);
        ParticleLeniaBucketThreaded particle_lenia(parameters, particles, 1e-4, 0.5, human_readable, 6);
        //ParticleLenia particle_lenia(parameters, particles);

        auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < 100; ++i) {
            particle_lenia.step();
        }
        auto end = std::chrono::steady_clock::now();
        long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        double time_per_step = ((double) duration) / 100;
        if (human_readable) {
            std::cout << "Total time: "
                      << duration
                      << "ms" << std::endl;
            std::cout << "Time per step: "
                      << time_per_step
                      << "ms" << std::endl;
        } else {
            std::cout << num_particles << '\t' << duration << '\t' << time_per_step << std::endl;
        }
        /*if (duration > max_time) {
            if (human_readable) {
                std::cout << "Execution time exceeded " << max_time << "ms, stopping..." << std::endl;
            }
            break;
        }*/
        if (num_particles > 100000) {
            num_particles += 100000;
        } else {
            num_particles *= 2;
        }
    }
}