#ifndef PARTICLE_LENIA_CPU_PARTICLELENIA_H
#define PARTICLE_LENIA_CPU_PARTICLELENIA_H

#include <vector>
#include "Parameters.h"
#include "Fields.h"

struct Point {
    double x;
    double y;

    double difference_norm(const Point &other) const;

    double norm() const;

    double normalize();

    Point operator+(const Point &other) const;

    Point operator-(const Point &other) const;

    Point &operator+=(const Point &other);

    Point &operator-=(const Point &other);

    Point operator*(double factor) const;

    Point &operator*=(double factor);
};

Point operator*(double factor, const Point &point);

class ParticleLenia {
public:
    explicit ParticleLenia(Parameters parameters, const std::vector<Point> &points);

    Parameters parameters;
    std::vector<Point> points;
    std::vector<Point> gradients;

    Fields calculate_fields(Point point) const;

    void step();

private:
    static double kernel_function(double x, double mu, double sigma2);

    Point gradient(const Point &point) const;
};


#endif //PARTICLE_LENIA_CPU_PARTICLELENIA_H
