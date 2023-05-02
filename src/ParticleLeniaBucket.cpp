//
// Created by adri on 31.01.23.
//

#include <complex>
#include <iostream>
#include "ParticleLeniaBucket.h"

ParticleLeniaBucket::ParticleLeniaBucket(Parameters parameters, std::vector<Point> points,
                                         double max_individual_step_error, double relative_bucket_size,
                                         bool debug_messages)
        : parameters(parameters) {
    number_particles = (int) points.size();
    max_distance = find_max_distance(max_individual_step_error, 1e-3, 1e3, debug_messages);
    max_distance_squared = std::pow(max_distance, 2);
    bucket_size = max_distance * relative_bucket_size;
    traversal_distance = std::ceil(max_distance / bucket_size);

    for (int i = 0; i < points.size(); ++i) {
        buckets[to_identifier(points[i])].points.push_back(
                {
                        i, points[i].x, points[i].y
                }
        );
        dU.push_back({0, 0});
        dR.push_back({0, 0});
        U.push_back(0);
    }
}

// in the following function we assume G'(U) = 1, since G'(U) <= 1
// as well as expecting distance > 1
double ParticleLeniaBucket::test_distance(double distance) const {
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

double ParticleLeniaBucket::find_max_distance(double max_error, double min_step_size, double max_step_size,
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

void ParticleLeniaBucket::gradient(const IdentifiablePoint &point) {
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

                        dR[point.id] += (dr * diff);
                        dR[other_point.id] -= (dr * diff);
                    }

                    double kernel = kernel_function(norm, parameters.mu_k, sigma_k2);
                    U[point.id] += kernel;
                    U[other_point.id] += kernel;
                    double du = (norm - parameters.mu_k) * kernel;
                    dU[point.id] += (du * diff);
                    dU[other_point.id] -= (du * diff);
                }
            }
        }
    }
}

void ParticleLeniaBucket::step() {
    double sigma_k2 = std::pow(parameters.sigma_k, 2);
    double sigma_g2 = std::pow(parameters.sigma_g, 2);

    std::fill(dU.begin(), dU.end(), Point{0, 0});
    std::fill(dR.begin(), dR.end(), Point{0, 0});
    std::fill(U.begin(), U.end(), kernel_function(0, parameters.mu_k, sigma_k2));

    std::map<BucketIdentifier, Bucket> updated_buckets;

    for (const auto &bucket: buckets) {
        for (IdentifiablePoint point: bucket.second.points) {
            gradient(point);
        }
    }

    for (const auto &bucket: buckets) {
        for (IdentifiablePoint point: bucket.second.points) {
            double u = U[point.id] * parameters.w_k;
            Point du = (-(2 * parameters.w_k / sigma_k2)) * dU[point.id];
            Point dr = dR[point.id];
            Point gradient =
            point -= parameters.dt * gradient;
            updated_buckets[to_identifier(point)].points.push_back(point);
        }
    }

    std::swap(buckets, updated_buckets);
}

double ParticleLeniaBucket::kernel_function(double x, double mu, double sigma2) {
    return std::exp(-std::pow(x - mu, 2) / sigma2);
}

ParticleLeniaBucket::BucketIdentifier ParticleLeniaBucket::to_identifier(const Point &point) const {
    return {
            std::ceil(point.x / bucket_size),
            std::ceil(point.y / bucket_size)
    };
}

ParticleLeniaBucket::BucketIdentifier ParticleLeniaBucket::to_identifier(const IdentifiablePoint &point) const {
    return {
            std::ceil(point.x / bucket_size),
            std::ceil(point.y / bucket_size)
    };
}


double IdentifiablePoint::difference_norm(const IdentifiablePoint &other) const {
    return std::sqrt(std::pow(x - other.x, 2) + std::pow(y - other.y, 2));
}

IdentifiablePoint &IdentifiablePoint::operator+=(const Point &other) {
    x += other.x;
    y += other.y;
    return *this;
}

IdentifiablePoint &IdentifiablePoint::operator-=(const Point &other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

double IdentifiablePoint::norm() const {
    return std::sqrt(std::pow(x, 2) + std::pow(y, 2));
}

double IdentifiablePoint::normalize() {
    double length = norm();
    x /= length;
    y /= length;
    return length;
}

IdentifiablePoint &IdentifiablePoint::operator*=(double factor) {
    x *= factor;
    y *= factor;
    return *this;
}

Point IdentifiablePoint::to_point() const {
    return {
            x,
            y
    };
}

Point IdentifiablePoint::operator-(const IdentifiablePoint &other) const {
    return {
            x - other.x,
            y - other.y
    };
}
