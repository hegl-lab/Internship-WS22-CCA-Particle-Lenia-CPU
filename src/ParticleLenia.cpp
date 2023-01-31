//
// Created by adri on 31.01.23.
//

#include <complex>
#include "ParticleLenia.h"

Fields ParticleLenia::calculate_fields(Point point) const {
    double u = 0;
    double r = 0;

    double sigma_k2 = std::pow(parameters.sigma_k, 2);

    for (auto &other_point: points) {
        double norm = point.difference_norm(other_point);
        u += kernel_function(norm, parameters.mu_k, sigma_k2);
        if (norm >= 1e-10) r += std::pow(std::max(1 - norm, 0.), 2);
    }

    r *= (parameters.c_rep * 0.5);
    u *= parameters.w_k;

    double g = kernel_function(u, parameters.mu_g, std::pow(parameters.sigma_g, 2));
    double e = r - g;

    return {
            u, r, g, e
    };
}


double ParticleLenia::kernel_function(double x, double mu, double sigma2) {
    return std::exp(-std::pow(x - mu, 2) / sigma2);
}

void ParticleLenia::step() {
    for (int i = 0; i < points.size(); ++i) {
        gradients[i] = gradient(points[i]);
    }
    for (int i = 0; i < points.size(); ++i) {
        points[i] -= parameters.dt * gradients[i];
    }
}

ParticleLenia::ParticleLenia(Parameters parameters, const std::vector<Point> &points) : parameters(parameters),
                                                                                        points(points) {
    for (int i = 0; i < points.size(); ++i) {
        gradients.push_back({0, 0});
    }
}

Point ParticleLenia::gradient(const Point &point) const {
    Point dU = {0, 0};
    Point dR = {0, 0};
    double sigma_k2 = std::pow(parameters.sigma_k, 2);
    double U = kernel_function(0, parameters.mu_k, sigma_k2);

    for (auto &other_point: points) {
        Point diff = other_point - point;
        double norm = diff.normalize();
        if (norm < 1e-5) continue;

        if (norm < 1) {
            double dr = -parameters.c_rep * (1 - norm);
            dR += dr * diff;
        }

        double kernel = kernel_function(norm, parameters.mu_k, sigma_k2);
        U += kernel;
        double du = (norm - parameters.mu_k) * kernel;
        dU += du * diff;
    }

    U *= parameters.w_k;
    dU *= -(2 * parameters.w_k / sigma_k2);

    double sigma_g2 = std::pow(parameters.sigma_g, 2);
    Point dG = (-2 * (U - parameters.mu_g) * kernel_function(U, parameters.mu_g, sigma_g2) / sigma_g2) * dU - dR;
    return dG;
}

double Point::difference_norm(const Point &other) const {
    return std::sqrt(std::pow(x - other.x, 2) + std::pow(y - other.y, 2));
}

Point Point::operator+(const Point &other) const {
    return {
            x + other.x,
            y + other.y
    };
}

Point Point::operator-(const Point &other) const {
    return {
            x - other.x,
            y - other.y
    };
}

Point &Point::operator+=(const Point &other) {
    x += other.x;
    y += other.y;
    return *this;
}

Point &Point::operator-=(const Point &other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Point Point::operator*(double factor) const {
    return {
            x * factor,
            y * factor
    };
}

double Point::norm() const {
    return std::sqrt(std::pow(x, 2) + std::pow(y, 2));
}

double Point::normalize() {
    double length = norm();
    x /= length;
    y /= length;
    return length;
}

Point &Point::operator*=(double factor) {
    x *= factor;
    y *= factor;
    return *this;
}

Point operator*(double factor, const Point &point) {
    return point * factor;
}
