#include <iostream>
#include <random>
#include <chrono>
#include "ParticleLenia.h"

int main() {
    int num_particles = 200;

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> distribution(0, 10);

    std::vector<Point> particles;
    for (int i = 0; i < num_particles; ++i)
        particles.push_back(Point{
                distribution(rng), distribution(rng)
        });

    Parameters parameters;
    ParticleLenia particle_lenia(parameters, particles);

    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < 2000; ++i) {
        particle_lenia.step();
    }
    auto end = std::chrono::steady_clock::now();
    long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double time_per_step = ((double) duration) / 2000;
    std::cout << "Total time: "
              << duration
              << "ms" << std::endl;
    std::cout << "Time per step: "
              << time_per_step
              << "ms" << std::endl;


    for (auto &point: particle_lenia.points) std::cout << point.x << '\t' << point.y << std::endl;

    return 0;
}
