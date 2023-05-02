#ifndef PARTICLE_LENIA_CPU_PARTICLELENIABUCKET_H
#define PARTICLE_LENIA_CPU_PARTICLELENIABUCKET_H


#include <map>
#include <vector>
#include "Parameters.h"
#include "ParticleLenia.h"

struct IdentifiablePoint {
    int id;
    double x;
    double y;

    double difference_norm(const IdentifiablePoint &other) const;

    double norm() const;

    double normalize();

    Point to_point() const;

    IdentifiablePoint &operator+=(const Point &other);

    IdentifiablePoint &operator-=(const Point &other);

    Point operator-(const IdentifiablePoint &other) const;

    IdentifiablePoint &operator*=(double factor);
};

class ParticleLeniaBucket {
public:
    using BucketIdentifier = std::pair<int, int>;
    struct Bucket {
        std::vector<IdentifiablePoint> points;
    };

    ParticleLeniaBucket(Parameters parameters, std::vector<Point> points, double max_individual_step_error,
                        double relative_bucket_size, bool debug_messages);

    void gradient(const IdentifiablePoint &point);

    void step();

    double max_distance;
    double max_distance_squared;
    double bucket_size;
    int traversal_distance;
    int number_particles;
    Parameters parameters;
    std::map<BucketIdentifier, Bucket> buckets;
    std::vector<Point> dU;
    std::vector<Point> dR;
    std::vector<double> U;
private:
    double test_distance(double distance) const;

    double find_max_distance(double max_error, double min_step_size, double max_step_size, bool debug_messages) const;

    BucketIdentifier to_identifier(const IdentifiablePoint &point) const;

    BucketIdentifier to_identifier(const Point &point) const;

    static double kernel_function(double x, double mu, double sigma2);
};


#endif //PARTICLE_LENIA_CPU_PARTICLELENIABUCKET_H
