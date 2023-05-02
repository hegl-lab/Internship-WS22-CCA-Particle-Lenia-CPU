#include <complex>
#include <iostream>
#include <thread>
#include "ParticleLeniaBucketThreaded.h"

ParticleLeniaBucketThreaded::ParticleLeniaBucketThreaded(Parameters parameters, std::vector<Point> points,
                                                         double max_individual_step_error, double relative_bucket_size,
                                                         bool debug_messages, int threads)
        : parameters(parameters), number_threads(threads) {
    max_distance = find_max_distance(max_individual_step_error, 1e-3, 1e3, debug_messages);
    max_distance_squared = std::pow(max_distance, 2);
    bucket_size = max_distance * relative_bucket_size;
    traversal_distance = std::ceil(max_distance / bucket_size);
    number_particles = (int) points.size();

    std::vector<Point> empty_points;
    std::vector<double> empty_double;

    for (int i = 0; i < points.size(); ++i) {
        buckets[to_identifier(points[i])].points.push_back(
                {
                        i, points[i].x, points[i].y
                }
        );
        this->points.push_back(
                {
                        i, points[i].x, points[i].y
                }
        );
        empty_points.push_back({0, 0});
        empty_double.push_back(0);
    }

    for (int i = 0; i < threads; ++i) {
        dU.push_back(empty_points);
        dR.push_back(empty_points);
        U.push_back(empty_double);
    }
}

// in the following function we assume G'(U) = 1, since G'(U) <= 1
// as well as expecting distance > 1
double ParticleLeniaBucketThreaded::test_distance(double distance) const {
    double sigma_k2 = std::pow(parameters.sigma_k, 2);
    double dU = (distance - parameters.mu_k) * kernel_function(distance, parameters.mu_k, sigma_k2);
    double U = kernel_function(0, parameters.mu_k, sigma_k2) +
               kernel_function(distance, parameters.mu_k, sigma_k2);

    U *= parameters.w_k;
    dU *= -(2 * parameters.w_k / sigma_k2);

    double sigma_g2 = std::pow(parameters.sigma_g, 2);
    double dg = std::abs(
            (-2 * ((number_particles - 1) * U - parameters.mu_g) * 1 / sigma_g2) * (number_particles - 1) * dU
    );
    return dg;
}

double ParticleLeniaBucketThreaded::find_max_distance(double max_error, double min_step_size, double max_step_size,
                                                      bool debug_messages) const {
    max_error /= parameters.dt;
    double size = 1.0;
    double step_size = max_step_size;
    int steps = 0;

    while (step_size >= min_step_size) {
        while (test_distance(size + step_size) > max_error) {
            size += step_size;
            ++steps;
        }
        step_size /= 2;
        ++steps;
    }
    size += step_size * 2;

    double error = test_distance(size + step_size);
    if (debug_messages)
        std::cout << "Found minimal bucket size of " << size << " with maximum step error of " << error * parameters.dt
                  << " in " << steps
                  << " steps." << std::endl;

    return size;
}

void ParticleLeniaBucketThreaded::gradient(const IdentifiablePoint &point, std::vector<Point> &dU_local,
                                           std::vector<Point> &dR_local, std::vector<double> &U_local) {
    double sigma_k2 = std::pow(parameters.sigma_k, 2);

    BucketIdentifier identifier = to_identifier(point);

    for (int i = -traversal_distance; i <= traversal_distance; ++i) {
        for (int j = -traversal_distance; j <= traversal_distance; ++j) {
            BucketIdentifier target_identifier{
                    identifier.first + i,
                    identifier.second + j
            };
            auto bucket_iterator = buckets.find(target_identifier);
            if (bucket_iterator != buckets.end()) {
                for (const auto &other_point: bucket_iterator->second.points) {
                    Point diff = other_point - point;
                    if (diff.x < 0) continue;
                    double norm = std::pow(diff.x, 2) + std::pow(diff.y, 2);
                    if (norm > max_distance_squared) continue;
                    norm = std::sqrt(norm);
                    if (norm < 1e-5) continue;
                    diff.x /= norm;
                    diff.y /= norm;

                    if (norm < 1) {
                        double dr = -parameters.c_rep * (1 - norm);

                        dR_local[point.id] += (dr * diff);
                        dR_local[other_point.id] -= (dr * diff);
                    }

                    double kernel = kernel_function(norm, parameters.mu_k, sigma_k2);
                    U_local[point.id] += kernel;
                    U_local[other_point.id] += kernel;
                    double du = (norm - parameters.mu_k) * kernel;
                    dU_local[point.id] += (du * diff);
                    dU_local[other_point.id] -= (du * diff);
                }
            }
        }
    }
}


void ParticleLeniaBucketThreaded::thread_step(int id, int start_index, int end_index) {
    double sigma_k2 = std::pow(parameters.sigma_k, 2);

    auto &dU_local = dU[id];
    auto &dR_local = dR[id];
    auto &U_local = U[id];

    std::fill(dU_local.begin(), dU_local.end(), Point{0, 0});
    std::fill(dR_local.begin(), dR_local.end(), Point{0, 0});
    std::fill(U_local.begin(), U_local.end(), kernel_function(0, parameters.mu_k, sigma_k2));

    for (int i = start_index; i < end_index; ++i) {
        gradient(points[i], dU_local, dR_local, U_local);
    }
}


void ParticleLeniaBucketThreaded::step() {
    int objects_per_thread = points.size() / number_threads;
    std::vector<std::thread> threads;
    for (int i = 0; i < number_threads; ++i) {
        threads.emplace_back(&ParticleLeniaBucketThreaded::thread_step, this, i, i * objects_per_thread,
                             ((i + 1) == number_threads ? points.size() : (i + 1) * objects_per_thread));
    }
    for (auto &thread: threads) thread.join();

    double sigma_k2 = std::pow(parameters.sigma_k, 2);
    double sigma_g2 = std::pow(parameters.sigma_g, 2);

    std::map<BucketIdentifier, Bucket> updated_buckets;
    std::vector<IdentifiablePoint> updated_points;

    for (const auto &bucket: buckets) {
        for (IdentifiablePoint point: bucket.second.points) {
            double u = aggregate_U(point.id) * parameters.w_k;
            Point du = (-(2 * parameters.w_k / sigma_k2)) * aggregate_dU(point.id);
            Point dr = aggregate_dR(point.id);
            Point gradient =
                    (-2 * (u - parameters.mu_g) * kernel_function(u, parameters.mu_g, sigma_g2) / sigma_g2) * du - dr;
            point -= parameters.dt * gradient;
            updated_buckets[to_identifier(point)].points.push_back(point);
            updated_points.push_back(point);
        }
    }

    std::swap(buckets, updated_buckets);
    std::swap(points, updated_points);
}

double ParticleLeniaBucketThreaded::kernel_function(double x, double mu, double sigma2) {
    return std::exp(-std::pow(x - mu, 2) / sigma2);
}

ParticleLeniaBucketThreaded::BucketIdentifier ParticleLeniaBucketThreaded::to_identifier(const Point &point) const {
    return {
            std::ceil(point.x / bucket_size),
            std::ceil(point.y / bucket_size)
    };
}

ParticleLeniaBucketThreaded::BucketIdentifier
ParticleLeniaBucketThreaded::to_identifier(const IdentifiablePoint &point) const {
    return {
            std::ceil(point.x / bucket_size),
            std::ceil(point.y / bucket_size)
    };
}

double ParticleLeniaBucketThreaded::aggregate_U(int id) const {
    double u = 0;
    for (auto &U_local: U) u += U_local[id];
    return u;
}

Point ParticleLeniaBucketThreaded::aggregate_dU(int id) const {
    Point du = {0, 0};
    for (auto &dU_local: dU) du += dU_local[id];
    return du;
}

Point ParticleLeniaBucketThreaded::aggregate_dR(int id) const {
    Point dr = {0, 0};
    for (auto &dR_local: dR) dr += dR_local[id];
    return dr;
}