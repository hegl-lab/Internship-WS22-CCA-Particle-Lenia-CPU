#include <random>
#include "raylib.h"
#include "ParticleLeniaBucket.h"
#include "ParticleLeniaBucketThreaded.h"

int main() {
    int num_particles = 200;

    double width = 1.1 * std::sqrt(num_particles);

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> distribution(0, width);

    double offset = -width * 0.25;
    width *= 1.5;

    std::vector<Point> particles;
    for (int i = 0; i < num_particles; ++i)
        particles.push_back(Point{
                distribution(rng), distribution(rng)
        });

    Parameters parameters;
    ParticleLeniaBucket particle_lenia(parameters, particles, 1e-10, 0.5, true);
    //ParticleLeniaBucketThreaded particle_lenia(parameters, particles, 1e-10, 0.5, true, 4);

    const int screenWidth = 800;
    const int screenHeight = 800;

    InitWindow(screenWidth, screenHeight, "Particle Lenia Bucket");

    SetTargetFPS(60);

    // Main render loop
    while (!WindowShouldClose()) {

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        ClearBackground(RAYWHITE);

        particle_lenia.step();
        for (const auto &bucket: particle_lenia.buckets) {
            for (IdentifiablePoint point: bucket.second.points) {
                Point position = point.to_point();
                position -= {offset, offset};
                position *= 1 / width;
                position *= screenWidth;
                DrawCircleV({(float) position.x, (float) position.y}, 5, BLACK);
            }
        }

        EndDrawing();
    }

    CloseWindow();

    return 0;
}
