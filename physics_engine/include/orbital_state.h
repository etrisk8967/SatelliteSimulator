// physics_engine/include/orbital_state.h
#pragma once

#include <string>
#include "vector_math.h"

namespace SatelliteSimulator {

    // Gravitational constant * Earth mass (m³/s²)
    constexpr double EARTH_MU = 3.986004418e14;
    
    // Earth radius in meters
    constexpr double EARTH_RADIUS = 6378137.0;

    // State representation using Cartesian coordinates
    class CartesianState {
    public:
        // Constructors
        CartesianState();
        CartesianState(const Vector3& position, const Vector3& velocity);
        
        // Getters & Setters
        Vector3 getPosition() const;
        Vector3 getVelocity() const;
        void setPosition(const Vector3& position);
        void setVelocity(const Vector3& velocity);
        
        // Utility functions
        double getKineticEnergy(double mass) const;
        double getPotentialEnergy(double mass) const;
        double getTotalEnergy(double mass) const;
        double getAngularMomentum(double mass) const;
        Vector3 getAngularMomentumVector() const;
        
        // String representation
        std::string toString() const;
        
    private:
        Vector3 position_; // in meters
        Vector3 velocity_; // in meters/second
    };
    
    // Keplerian Orbital Elements
    class OrbitalElements {
    public:
        // Constructors
        OrbitalElements();
        OrbitalElements(double semiMajorAxis, double eccentricity, double inclination,
                        double raan, double argOfPerigee, double trueAnomaly);
        
        // Getters & Setters
        double getSemiMajorAxis() const;
        double getEccentricity() const;
        double getInclination() const;
        double getRightAscension() const;
        double getArgumentOfPerigee() const;
        double getTrueAnomaly() const;
        
        void setSemiMajorAxis(double semiMajorAxis);
        void setEccentricity(double eccentricity);
        void setInclination(double inclination);
        void setRightAscension(double raan);
        void setArgumentOfPerigee(double argOfPerigee);
        void setTrueAnomaly(double trueAnomaly);
        
        // Derived properties
        double getPeriapsis() const;
        double getApoapsis() const;
        double getSemiLatusRectum() const;
        double getPeriod() const;  // Orbital period in seconds
        
        // String representation
        std::string toString() const;
        
    private:
        double semiMajorAxis_;   // a (meters)
        double eccentricity_;    // e (dimensionless)
        double inclination_;     // i (radians)
        double raan_;            // Ω (radians) - Right Ascension of the Ascending Node
        double argOfPerigee_;    // ω (radians) - Argument of Perigee
        double trueAnomaly_;     // ν (radians) - True Anomaly
    };
    
    // Conversion functions between state representations
    CartesianState orbitalElementsToCartesian(const OrbitalElements& elements, double mu = EARTH_MU);
    OrbitalElements cartesianToOrbitalElements(const CartesianState& state, double mu = EARTH_MU);
    
    // Anomaly conversions
    double meanToEccentricAnomaly(double meanAnomaly, double eccentricity, double tolerance = 1e-10);
    double eccentricToTrueAnomaly(double eccentricAnomaly, double eccentricity);
    double trueToEccentricAnomaly(double trueAnomaly, double eccentricity);
    double eccentricToMeanAnomaly(double eccentricAnomaly, double eccentricity);
    double trueToMeanAnomaly(double trueAnomaly, double eccentricity);
    double meanToTrueAnomaly(double meanAnomaly, double eccentricity, double tolerance = 1e-10);
}