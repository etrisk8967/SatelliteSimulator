// physics_engine/include/vector_math.h
#pragma once

#include <Eigen/Dense>
#include <string>

namespace SatelliteSimulator {
    // Type aliases for clarity
    using Vector3 = Eigen::Vector3d;
    using Matrix3 = Eigen::Matrix3d;

    // Basic vector operations
    Vector3 addVectors(const Vector3& v1, const Vector3& v2);
    Vector3 subtractVectors(const Vector3& v1, const Vector3& v2);
    Vector3 scaleVector(const Vector3& v, double scalar);
    double dotProduct(const Vector3& v1, const Vector3& v2);
    Vector3 crossProduct(const Vector3& v1, const Vector3& v2);

    // Vector properties
    double magnitude(const Vector3& v);
    Vector3 normalize(const Vector3& v);
    double angleBetweenVectors(const Vector3& v1, const Vector3& v2);

    // Coordinate transformations
    enum class CoordinateSystem {
        CARTESIAN,
        SPHERICAL,
        CYLINDRICAL
    };

    Vector3 transformCoordinates(const Vector3& v, 
                                CoordinateSystem from, 
                                CoordinateSystem to);

    // Helper for spherical coordinates (r, theta, phi)
    // theta is the polar angle (from z-axis), phi is the azimuthal angle (in x-y plane)
    Vector3 cartesianToSpherical(const Vector3& cartesian);
    Vector3 sphericalToCartesian(const Vector3& spherical);

    // Helper for cylindrical coordinates (rho, phi, z)
    Vector3 cartesianToCylindrical(const Vector3& cartesian);
    Vector3 cylindricalToCartesian(const Vector3& cylindrical);

    // Rotation functions
    Matrix3 rotationMatrixX(double angle);  // Angle in radians
    Matrix3 rotationMatrixY(double angle);
    Matrix3 rotationMatrixZ(double angle);
    Vector3 rotateVector(const Vector3& v, const Matrix3& rotationMatrix);
    
    // String formatting for debugging
    std::string vectorToString(const Vector3& v);
    std::string matrixToString(const Matrix3& m);
}