// physics_engine/src/vector_math.cpp
#include "../include/vector_math.h"
#include <cmath>
#include <sstream>
#include <iomanip>

namespace SatelliteSimulator {

    Vector3 addVectors(const Vector3& v1, const Vector3& v2) {
        return v1 + v2;
    }

    Vector3 subtractVectors(const Vector3& v1, const Vector3& v2) {
        return v1 - v2;
    }

    Vector3 scaleVector(const Vector3& v, double scalar) {
        return v * scalar;
    }

    double dotProduct(const Vector3& v1, const Vector3& v2) {
        return v1.dot(v2);
    }

    Vector3 crossProduct(const Vector3& v1, const Vector3& v2) {
        return v1.cross(v2);
    }

    double magnitude(const Vector3& v) {
        return v.norm();
    }

    Vector3 normalize(const Vector3& v) {
        double mag = magnitude(v);
        if (mag < 1e-10) {  // Prevent division by zero
            return Vector3::Zero();
        }
        return v / mag;
    }

    double angleBetweenVectors(const Vector3& v1, const Vector3& v2) {
        double dot = dotProduct(v1, v2);
        double mag1 = magnitude(v1);
        double mag2 = magnitude(v2);
        
        if (mag1 < 1e-10 || mag2 < 1e-10) {
            return 0.0;  // Undefined angle if either vector is zero
        }
        
        // Clamp to prevent numerical errors
        double cosAngle = std::clamp(dot / (mag1 * mag2), -1.0, 1.0);
        return std::acos(cosAngle);
    }

    Vector3 cartesianToSpherical(const Vector3& cartesian) {
        double x = cartesian.x();
        double y = cartesian.y();
        double z = cartesian.z();
        
        double r = std::sqrt(x*x + y*y + z*z);
        double theta = 0.0;
        double phi = 0.0;
        
        if (r > 1e-10) {
            theta = std::acos(z / r);  // polar angle (from z-axis)
            phi = std::atan2(y, x);    // azimuthal angle (in x-y plane)
        }
        
        return Vector3(r, theta, phi);
    }

    Vector3 sphericalToCartesian(const Vector3& spherical) {
        double r = spherical.x();
        double theta = spherical.y();  // polar angle
        double phi = spherical.z();    // azimuthal angle
        
        double x = r * std::sin(theta) * std::cos(phi);
        double y = r * std::sin(theta) * std::sin(phi);
        double z = r * std::cos(theta);
        
        return Vector3(x, y, z);
    }

    Vector3 cartesianToCylindrical(const Vector3& cartesian) {
        double x = cartesian.x();
        double y = cartesian.y();
        double z = cartesian.z();
        
        double rho = std::sqrt(x*x + y*y);
        double phi = std::atan2(y, x);
        
        return Vector3(rho, phi, z);
    }

    Vector3 cylindricalToCartesian(const Vector3& cylindrical) {
        double rho = cylindrical.x();
        double phi = cylindrical.y();
        double z = cylindrical.z();
        
        double x = rho * std::cos(phi);
        double y = rho * std::sin(phi);
        
        return Vector3(x, y, z);
    }

    Vector3 transformCoordinates(const Vector3& v, 
                                CoordinateSystem from, 
                                CoordinateSystem to) {
        // First convert to Cartesian if not already
        Vector3 cartesian;
        switch(from) {
            case CoordinateSystem::CARTESIAN:
                cartesian = v;
                break;
            case CoordinateSystem::SPHERICAL:
                cartesian = sphericalToCartesian(v);
                break;
            case CoordinateSystem::CYLINDRICAL:
                cartesian = cylindricalToCartesian(v);
                break;
        }
        
        // Then convert from Cartesian to target system
        switch(to) {
            case CoordinateSystem::CARTESIAN:
                return cartesian;
            case CoordinateSystem::SPHERICAL:
                return cartesianToSpherical(cartesian);
            case CoordinateSystem::CYLINDRICAL:
                return cartesianToCylindrical(cartesian);
        }
        
        // This should never happen, but compiler might warn without it
        return Vector3::Zero();
    }

    Matrix3 rotationMatrixX(double angle) {
        Matrix3 rot;
        rot << 1, 0, 0,
               0, std::cos(angle), -std::sin(angle),
               0, std::sin(angle), std::cos(angle);
        return rot;
    }

    Matrix3 rotationMatrixY(double angle) {
        Matrix3 rot;
        rot << std::cos(angle), 0, std::sin(angle),
               0, 1, 0,
               -std::sin(angle), 0, std::cos(angle);
        return rot;
    }

    Matrix3 rotationMatrixZ(double angle) {
        Matrix3 rot;
        rot << std::cos(angle), -std::sin(angle), 0,
               std::sin(angle), std::cos(angle), 0,
               0, 0, 1;
        return rot;
    }

    Vector3 rotateVector(const Vector3& v, const Matrix3& rotationMatrix) {
        return rotationMatrix * v;
    }

    std::string vectorToString(const Vector3& v) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
        return oss.str();
    }

    std::string matrixToString(const Matrix3& m) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << std::endl;
        for (int i = 0; i < 3; ++i) {
            oss << "[ ";
            for (int j = 0; j < 3; ++j) {
                oss << m(i, j) << " ";
            }
            oss << "]" << std::endl;
        }
        return oss.str();
    }
}