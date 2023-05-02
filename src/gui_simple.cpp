#include <random>
#include "raylib.h"
#include "ParticleLeniaBucket.h"

int main() {
    int num_particles = 2000;

    double width = 0.7 * std::sqrt(num_particles);

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
    ParticleLenia particle_lenia(parameters, particles);

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
        for (Point position: particle_lenia.points) {
            position -= {offset, offset};
            position *= 1 / width;
            position *= screenWidth;
            DrawCircleV({(float) position.x, (float) position.y}, 5, BLACK);
        }

        EndDrawing();
    }

    CloseWindow();

    return 0;
}