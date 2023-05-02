#ifndef PARTICLE_LENIA_CPU_PARTICLELENIABUCKETTHREADED_H
#define PARTICLE_LENIA_CPU_PARTICLELENIABUCKETTHREADED_H


#include <map>
#include <vector>
#include "Parameters.h"
#include "ParticleLenia.h"
#include "ParticleLeniaBucket.h"

class ParticleLeniaBucketThreaded {
public:
    using BucketIdentifier = std::pair<int, int>;
    struct Bucket {
        std::vector<IdentifiablePoint> points;
    };

    ParticleLeniaBucketThreaded(Parameters parameters, std::vector<Point> points, double max_individual_step_error,
                                double relative_bucket_size, bool debug_messages, int threads);

    void gradient(const IdentifiablePoint &point, std::vector<Point> &dU_local,
                  std::vector<Point> &dR_local, std::vector<double> &U_local);

    void step();

    double max_distance;
    double max_distance_squared;
    double bucket_size;
    int traversal_distance;
    int number_particles;
    int number_threads;
    Parameters parameters;
    std::map<BucketIdentifier, Bucket> buckets;
    std::vector<IdentifiablePoint> points;
    std::vector<std::vector<Point>> dU;
    std::vector<std::vector<Point>> dR;
    std::vector<std::vector<double>> U;
private:
    double test_distance(double distance) const;

    double find_max_distance(double max_error, double min_step_size, double max_step_size, bool debug_messages) const;

    BucketIdentifier to_identifier(const IdentifiablePoint &point) const;

    BucketIdentifier to_identifier(const Point &point) const;

    void thread_step(int id, int start_index, int end_index);

    double aggregate_U(int id) const;
    Point aggregate_dU(int id) const;
    Point aggregate_dR(int id) const;

    static double kernel_function(double x, double mu, double sigma2);
};


#endif //PARTICLE_LENIA_CPU_PARTICLELENIABUCKETTHREADED_H
