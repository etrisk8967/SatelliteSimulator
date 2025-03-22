// physics_engine/src/orbital_state.cpp
#define _USE_MATH_DEFINES  // This must come before including cmath
#include <cmath>
#include "../include/orbital_state.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace SatelliteSimulator {

    //--------------------------------------------------------
    // CartesianState Implementation
    //--------------------------------------------------------
    
    CartesianState::CartesianState() 
        : position_(Vector3::Zero()), velocity_(Vector3::Zero()) {
    }
    
    CartesianState::CartesianState(const Vector3& position, const Vector3& velocity) 
        : position_(position), velocity_(velocity) {
    }
    
    Vector3 CartesianState::getPosition() const {
        return position_;
    }
    
    Vector3 CartesianState::getVelocity() const {
        return velocity_;
    }
    
    void CartesianState::setPosition(const Vector3& position) {
        position_ = position;
    }
    
    void CartesianState::setVelocity(const Vector3& velocity) {
        velocity_ = velocity;
    }
    
    double CartesianState::getKineticEnergy(double mass) const {
        return 0.5 * mass * velocity_.squaredNorm();
    }
    
    double CartesianState::getPotentialEnergy(double mass) const {
        double r = position_.norm();
        if (r < 1e-10) {
            throw std::runtime_error("Position too close to origin, potential energy undefined");
        }
        return -EARTH_MU * mass / r;
    }
    
    double CartesianState::getTotalEnergy(double mass) const {
        return getKineticEnergy(mass) + getPotentialEnergy(mass);
    }
    
    Vector3 CartesianState::getAngularMomentumVector() const {
        return position_.cross(velocity_);
    }
    
    double CartesianState::getAngularMomentum(double mass) const {
        return mass * getAngularMomentumVector().norm();
    }
    
    std::string CartesianState::toString() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2);
        oss << "Position: " << vectorToString(position_) << " m" << std::endl;
        oss << "Velocity: " << vectorToString(velocity_) << " m/s";
        return oss.str();
    }
    
    //--------------------------------------------------------
    // OrbitalElements Implementation
    //--------------------------------------------------------
    
    OrbitalElements::OrbitalElements() 
        : semiMajorAxis_(0.0), eccentricity_(0.0), inclination_(0.0),
          raan_(0.0), argOfPerigee_(0.0), trueAnomaly_(0.0) {
    }
    
    OrbitalElements::OrbitalElements(double semiMajorAxis, double eccentricity, 
                                     double inclination, double raan, 
                                     double argOfPerigee, double trueAnomaly) 
        : semiMajorAxis_(semiMajorAxis), eccentricity_(eccentricity), 
          inclination_(inclination), raan_(raan), 
          argOfPerigee_(argOfPerigee), trueAnomaly_(trueAnomaly) {
          
        // Validate inputs
        if (semiMajorAxis <= 0 && eccentricity < 1.0) {
            throw std::invalid_argument("Semi-major axis must be positive for elliptical orbits");
        }
        
        if (eccentricity < 0) {
            throw std::invalid_argument("Eccentricity cannot be negative");
        }
        
        // For hyperbolic orbits, semi-major axis should be negative
        if (eccentricity > 1.0 && semiMajorAxis > 0) {
            semiMajorAxis_ = -semiMajorAxis;
        }
    }
    
    double OrbitalElements::getSemiMajorAxis() const {
        return semiMajorAxis_;
    }
    
    double OrbitalElements::getEccentricity() const {
        return eccentricity_;
    }
    
    double OrbitalElements::getInclination() const {
        return inclination_;
    }
    
    double OrbitalElements::getRightAscension() const {
        return raan_;
    }
    
    double OrbitalElements::getArgumentOfPerigee() const {
        return argOfPerigee_;
    }
    
    double OrbitalElements::getTrueAnomaly() const {
        return trueAnomaly_;
    }
    
    void OrbitalElements::setSemiMajorAxis(double semiMajorAxis) {
        if (semiMajorAxis <= 0 && eccentricity_ < 1.0) {
            throw std::invalid_argument("Semi-major axis must be positive for elliptical orbits");
        }
        
        // For hyperbolic orbits, semi-major axis should be negative
        if (eccentricity_ > 1.0 && semiMajorAxis > 0) {
            semiMajorAxis_ = -semiMajorAxis;
        } else {
            semiMajorAxis_ = semiMajorAxis;
        }
    }
    
    void OrbitalElements::setEccentricity(double eccentricity) {
        if (eccentricity < 0) {
            throw std::invalid_argument("Eccentricity cannot be negative");
        }
        eccentricity_ = eccentricity;
        
        // If orbit type changes (elliptical <-> hyperbolic), adjust semi-major axis sign
        if (eccentricity > 1.0 && semiMajorAxis_ > 0) {
            semiMajorAxis_ = -semiMajorAxis_;
        } else if (eccentricity < 1.0 && semiMajorAxis_ < 0) {
            semiMajorAxis_ = -semiMajorAxis_;
        }
    }
    
    void OrbitalElements::setInclination(double inclination) {
        inclination_ = inclination;
    }
    
    void OrbitalElements::setRightAscension(double raan) {
        raan_ = raan;
    }
    
    void OrbitalElements::setArgumentOfPerigee(double argOfPerigee) {
        argOfPerigee_ = argOfPerigee;
    }
    
    void OrbitalElements::setTrueAnomaly(double trueAnomaly) {
        trueAnomaly_ = trueAnomaly;
    }
    
    double OrbitalElements::getSemiLatusRectum() const {
        return semiMajorAxis_ * (1.0 - eccentricity_ * eccentricity_);
    }
    
    double OrbitalElements::getPeriapsis() const {
        return semiMajorAxis_ * (1.0 - eccentricity_);
    }
    
    double OrbitalElements::getApoapsis() const {
        if (eccentricity_ >= 1.0) {
            throw std::runtime_error("Apoapsis undefined for parabolic and hyperbolic orbits");
        }
        return semiMajorAxis_ * (1.0 + eccentricity_);
    }
    
    double OrbitalElements::getPeriod() const {
        if (eccentricity_ >= 1.0) {
            throw std::runtime_error("Period undefined for parabolic and hyperbolic orbits");
        }
        return 2.0 * M_PI * std::sqrt(std::pow(semiMajorAxis_, 3) / EARTH_MU);
    }
    
    std::string OrbitalElements::toString() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2);
        oss << "Semi-major axis: " << semiMajorAxis_ << " m" << std::endl;
        oss << "Eccentricity: " << eccentricity_ << std::endl;
        oss << "Inclination: " << inclination_ * 180.0 / M_PI << " deg" << std::endl;
        oss << "RAAN: " << raan_ * 180.0 / M_PI << " deg" << std::endl;
        oss << "Argument of Perigee: " << argOfPerigee_ * 180.0 / M_PI << " deg" << std::endl;
        oss << "True Anomaly: " << trueAnomaly_ * 180.0 / M_PI << " deg";
        
        if (eccentricity_ < 1.0) {
            oss << std::endl << "Periapsis: " << getPeriapsis() << " m" << std::endl;
            oss << "Apoapsis: " << getApoapsis() << " m" << std::endl;
            oss << "Period: " << getPeriod() / 60.0 << " min";
        } else {
            oss << std::endl << "Periapsis: " << getPeriapsis() << " m";
        }
        
        return oss.str();
    }
    
    //--------------------------------------------------------
    // Anomaly Conversion Functions
    //--------------------------------------------------------
    
    double meanToEccentricAnomaly(double meanAnomaly, double eccentricity, double tolerance) {
        // Normalize mean anomaly to [0, 2π]
        meanAnomaly = std::fmod(meanAnomaly, 2.0 * M_PI);
        if (meanAnomaly < 0) {
            meanAnomaly += 2.0 * M_PI;
        }
        
        if (eccentricity < 1.0) {  // Elliptical orbits
            // Initial guess (for small eccentricities, M ≈ E)
            double E = meanAnomaly;
            
            // Newton-Raphson iteration to solve Kepler's equation: M = E - e*sin(E)
            double delta = 1.0;
            int iterations = 0;
            
            while (std::abs(delta) > tolerance && iterations < 100) {
                delta = (E - eccentricity * std::sin(E) - meanAnomaly) / (1.0 - eccentricity * std::cos(E));
                E -= delta;
                iterations++;
            }
            
            if (iterations >= 100) {
                throw std::runtime_error("Failed to converge in meanToEccentricAnomaly");
            }
            
            return E;
        } else if (eccentricity > 1.0) {  // Hyperbolic orbits
            // Initial guess
            double F = std::log(2.0 * meanAnomaly / eccentricity + 1.8);
            
            // Newton-Raphson iteration to solve hyperbolic Kepler's equation: M = e*sinh(F) - F
            double delta = 1.0;
            int iterations = 0;
            
            while (std::abs(delta) > tolerance && iterations < 100) {
                delta = (eccentricity * std::sinh(F) - F - meanAnomaly) / (eccentricity * std::cosh(F) - 1.0);
                F -= delta;
                iterations++;
            }
            
            if (iterations >= 100) {
                throw std::runtime_error("Failed to converge in meanToEccentricAnomaly for hyperbolic orbit");
            }
            
            return F;
        } else {  // Parabolic orbits (e = 1)
            throw std::invalid_argument("Parabolic orbits (e=1) not supported in meanToEccentricAnomaly");
        }
    }
    
    double eccentricToTrueAnomaly(double eccentricAnomaly, double eccentricity) {
        if (eccentricity < 1.0) {  // Elliptical orbits
            double cosE = std::cos(eccentricAnomaly);
            double sinE = std::sin(eccentricAnomaly);
            
            // Calculate true anomaly using the formula: tan(v/2) = sqrt((1+e)/(1-e)) * tan(E/2)
            double tanHalfV = std::sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) * std::tan(eccentricAnomaly / 2.0);
            return 2.0 * std::atan(tanHalfV);
        } else if (eccentricity > 1.0) {  // Hyperbolic orbits
            // F is the hyperbolic eccentric anomaly
            double F = eccentricAnomaly;
            
            // Calculate true anomaly using the formula: tan(v/2) = sqrt((e+1)/(e-1)) * tanh(F/2)
            double tanHalfV = std::sqrt((eccentricity + 1.0) / (eccentricity - 1.0)) * std::tanh(F / 2.0);
            return 2.0 * std::atan(tanHalfV);
        } else {  // Parabolic orbits (e = 1)
            throw std::invalid_argument("Parabolic orbits (e=1) not supported in eccentricToTrueAnomaly");
        }
    }
    
    double trueToEccentricAnomaly(double trueAnomaly, double eccentricity) {
        if (eccentricity < 1.0) {  // Elliptical orbits
            // Calculate eccentric anomaly using the formula: tan(E/2) = sqrt((1-e)/(1+e)) * tan(v/2)
            double tanHalfE = std::sqrt((1.0 - eccentricity) / (1.0 + eccentricity)) * std::tan(trueAnomaly / 2.0);
            return 2.0 * std::atan(tanHalfE);
        } else if (eccentricity > 1.0) {  // Hyperbolic orbits
            // Calculate hyperbolic eccentric anomaly
            double tanHalfF = std::sqrt((eccentricity - 1.0) / (eccentricity + 1.0)) * std::tan(trueAnomaly / 2.0);
            return 2.0 * std::atanh(tanHalfF);
        } else {  // Parabolic orbits (e = 1)
            throw std::invalid_argument("Parabolic orbits (e=1) not supported in trueToEccentricAnomaly");
        }
    }
    
    double eccentricToMeanAnomaly(double eccentricAnomaly, double eccentricity) {
        if (eccentricity < 1.0) {  // Elliptical orbits
            return eccentricAnomaly - eccentricity * std::sin(eccentricAnomaly);
        } else if (eccentricity > 1.0) {  // Hyperbolic orbits
            // F is the hyperbolic eccentric anomaly
            double F = eccentricAnomaly;
            return eccentricity * std::sinh(F) - F;
        } else {  // Parabolic orbits (e = 1)
            throw std::invalid_argument("Parabolic orbits (e=1) not supported in eccentricToMeanAnomaly");
        }
    }
    
    double trueToMeanAnomaly(double trueAnomaly, double eccentricity) {
        double eccentricAnomaly = trueToEccentricAnomaly(trueAnomaly, eccentricity);
        return eccentricToMeanAnomaly(eccentricAnomaly, eccentricity);
    }
    
    double meanToTrueAnomaly(double meanAnomaly, double eccentricity, double tolerance) {
        double eccentricAnomaly = meanToEccentricAnomaly(meanAnomaly, eccentricity, tolerance);
        return eccentricToTrueAnomaly(eccentricAnomaly, eccentricity);
    }
    
    //--------------------------------------------------------
    // Conversion Functions
    //--------------------------------------------------------
    
    CartesianState orbitalElementsToCartesian(const OrbitalElements& elements, double mu) {
        // Extract orbital elements
        double a = elements.getSemiMajorAxis();
        double e = elements.getEccentricity();
        double i = elements.getInclination();
        double Omega = elements.getRightAscension();
        double omega = elements.getArgumentOfPerigee();
        double nu = elements.getTrueAnomaly();
        
        // Calculate the semi-latus rectum
        double p = a * (1.0 - e * e);
        
        // Handle hyperbolic orbits
        if (e > 1.0) {
            p = -a * (e * e - 1.0);
        }
        
        // Calculate distance from focus
        double r = p / (1.0 + e * std::cos(nu));
        
        // Calculate position in orbital plane (perifocal coordinates)
        Vector3 position_orbital;
        position_orbital.x() = r * std::cos(nu);
        position_orbital.y() = r * std::sin(nu);
        position_orbital.z() = 0.0;
        
        // Calculate velocity in orbital plane
        double sqrt_mu_p = std::sqrt(mu / p);
        Vector3 velocity_orbital;
        velocity_orbital.x() = -sqrt_mu_p * std::sin(nu);
        velocity_orbital.y() = sqrt_mu_p * (e + std::cos(nu));
        velocity_orbital.z() = 0.0;
        
        // Create rotation matrices for orbital to inertial transformation
        Matrix3 R3_Omega = rotationMatrixZ(Omega);
        Matrix3 R1_i = rotationMatrixX(i);
        Matrix3 R3_omega = rotationMatrixZ(omega);
        
        // Combined rotation matrix (applying rotations in the correct order)
        Matrix3 rotMatrix = R3_Omega * R1_i * R3_omega;
        
        // Transform vectors from orbital to inertial frame
        Vector3 position = rotMatrix * position_orbital;
        Vector3 velocity = rotMatrix * velocity_orbital;
        
        return CartesianState(position, velocity);
    }
    
    OrbitalElements cartesianToOrbitalElements(const CartesianState& state, double mu) {
        // Extract position and velocity vectors
        Vector3 r = state.getPosition();
        Vector3 v = state.getVelocity();
        
        // Calculate the magnitude of position and velocity
        double rmag = r.norm();
        double vmag = v.norm();
        
        // Angular momentum vector
        Vector3 h = r.cross(v);
        double hmag = h.norm();
        
        // Vector pointing towards ascending node (n)
        Vector3 z(0, 0, 1);
        Vector3 n = z.cross(h);
        double nmag = n.norm();
        
        // Eccentricity vector
        Vector3 e_vec = ((vmag * vmag - mu / rmag) * r - r.dot(v) * v) / mu;
        double e = e_vec.norm();
        
        // Semi-major axis
        double energy = vmag * vmag / 2.0 - mu / rmag;
        double a;
        
        if (std::abs(energy) < 1e-10) {
            // Parabolic orbit (e ≈ 1)
            a = std::numeric_limits<double>::infinity();
            e = 1.0;
        } else if (energy < 0) {
            // Elliptical orbit (e < 1)
            a = -mu / (2.0 * energy);
        } else {
            // Hyperbolic orbit (e > 1)
            a = mu / (2.0 * energy);
            a = -a;  // Negative semi-major axis for hyperbolic orbits
        }
        
        // Inclination
        double i = std::acos(h.z() / hmag);
        
        // Right Ascension of the Ascending Node (RAAN)
        double Omega;
        if (nmag < 1e-10) {
            // For equatorial orbits, RAAN is undefined, set to zero
            Omega = 0.0;
        } else {
            Omega = std::acos(n.x() / nmag);
            if (n.y() < 0) {
                Omega = 2.0 * M_PI - Omega;
            }
        }
        
        // Argument of Periapsis
        double omega;
        if (nmag < 1e-10) {
            // For equatorial orbits, use x-axis as reference
            omega = std::acos(e_vec.x() / e);
            if (e_vec.y() < 0) {
                omega = 2.0 * M_PI - omega;
            }
        } else if (e < 1e-10) {
            // For circular orbits, argument of periapsis is undefined, set to zero
            omega = 0.0;
        } else {
            omega = std::acos(n.dot(e_vec) / (nmag * e));
            if (e_vec.z() < 0) {
                omega = 2.0 * M_PI - omega;
            }
        }
        
        // True Anomaly
        double nu;
        if (e < 1e-10) {
            // For circular orbits, measure from ascending node
            if (nmag < 1e-10) {
                // For circular equatorial orbits, measure from x-axis
                nu = std::acos(r.x() / rmag);
                if (r.y() < 0) {
                    nu = 2.0 * M_PI - nu;
                }
            } else {
                nu = std::acos(n.dot(r) / (nmag * rmag));
                if (r.dot(n.cross(h)) < 0) {
                    nu = 2.0 * M_PI - nu;
                }
            }
        } else {
            nu = std::acos(e_vec.dot(r) / (e * rmag));
            if (r.dot(v) < 0) {
                nu = 2.0 * M_PI - nu;
            }
        }
        
        return OrbitalElements(a, e, i, Omega, omega, nu);
    }
}