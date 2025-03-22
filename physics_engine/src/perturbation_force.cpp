// physics_engine/src/perturbation_force.cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdexcept>
#include "../include/perturbation_force.h"

namespace SatelliteSimulator {

    //--------------------------------------------------------
    // J2 Perturbation Implementation
    //--------------------------------------------------------
    
    J2Perturbation::J2Perturbation(double j2Value)
        : j2Value_(j2Value) {
    }
    
    Vector3 J2Perturbation::calculateAcceleration(const CartesianState& state, double time) const {
        if (!enabled_) {
            return Vector3::Zero();
        }
        
        // Get position vector
        Vector3 position = state.getPosition();
        double x = position.x();
        double y = position.y();
        double z = position.z();
        
        // Calculate squared distance
        double r2 = position.squaredNorm();
        double r = std::sqrt(r2);
        
        if (r < 1e-10) {
            return Vector3::Zero(); // Avoid division by zero
        }
        
        // Calculate multiplier for J2 perturbation
        double j2Multiplier = 1.5 * j2Value_ * EARTH_MU * std::pow(EARTH_RADIUS, 2) / std::pow(r, 4);
        
        // Calculate acceleration components
        // See formula: a_J2 = (3J2·μ·R²)/(2r⁵) · [x(5z²/r² - 1), y(5z²/r² - 1), z(5z²/r² - 3)]
        double zr2 = z * z / r2;
        
        double ax = j2Multiplier * x * (5.0 * zr2 - 1.0);
        double ay = j2Multiplier * y * (5.0 * zr2 - 1.0);
        double az = j2Multiplier * z * (5.0 * zr2 - 3.0);
        
        return Vector3(ax, ay, az);
    }
    
    //--------------------------------------------------------
    // Atmospheric Drag Implementation
    //--------------------------------------------------------
    
    AtmosphericDrag::AtmosphericDrag(double dragCoefficient, double crossSectionalArea, double satelliteMass)
        : dragCoefficient_(dragCoefficient),
          crossSectionalArea_(crossSectionalArea),
          satelliteMass_(satelliteMass) {
        
        if (satelliteMass_ <= 0.0) {
            throw std::invalid_argument("Satellite mass must be positive");
        }
        
        if (crossSectionalArea_ <= 0.0) {
            throw std::invalid_argument("Cross-sectional area must be positive");
        }
        
        if (dragCoefficient_ < 0.0) {
            throw std::invalid_argument("Drag coefficient cannot be negative");
        }
    }
    
    Vector3 AtmosphericDrag::calculateAcceleration(const CartesianState& state, double time) const {
        if (!enabled_) {
            return Vector3::Zero();
        }
        
        // Get position and velocity
        Vector3 position = state.getPosition();
        Vector3 velocity = state.getVelocity();
        
        // Calculate altitude
        double altitude = position.norm() - EARTH_RADIUS;
        
        // No drag above certain altitude (e.g., 2000 km)
        if (altitude > 2000000.0) {
            return Vector3::Zero();
        }
        
        // Calculate atmospheric density
        double density = calculateAtmosphericDensity(altitude);
        
        // Calculate velocity magnitude
        double velocityMag = velocity.norm();
        
        if (velocityMag < 1e-10) {
            return Vector3::Zero(); // No drag for negligible velocity
        }
        
        // Calculate drag force: F = -0.5 * ρ * Cd * A * v² * v̂
        // Acceleration: a = F/m
        double dragFactor = -0.5 * density * dragCoefficient_ * crossSectionalArea_ * velocityMag / satelliteMass_;
        
        // Return acceleration vector (in the opposite direction of velocity)
        return dragFactor * velocity;
    }
    
    double AtmosphericDrag::calculateAtmosphericDensity(double altitude) const {
        // Simple exponential atmospheric model
        // ρ(h) = ρ₀ * exp(-h/H)
        // ρ₀ = density at sea level (kg/m³)
        // H = scale height (m)
        
        // Constants for Earth's atmosphere
        double rho0 = 1.225;      // Sea level density (kg/m³)
        double scaleHeight = 8500.0;  // Scale height (m)
        
        // Different density models for different altitude regions
        if (altitude < 0.0) {
            // Below sea level (shouldn't happen in normal simulations)
            return rho0;
        } else if (altitude <= 100000.0) {
            // 0-100 km: Exponential decay
            return rho0 * std::exp(-altitude / scaleHeight);
        } else if (altitude <= 500000.0) {
            // 100-500 km: Different scale height
            double rho100 = rho0 * std::exp(-100000.0 / scaleHeight);
            double scaleHeight2 = 60000.0;
            return rho100 * std::exp(-(altitude - 100000.0) / scaleHeight2);
        } else {
            // Above 500 km: Very low constant density
            return 1e-15;
        }
    }
    
    //--------------------------------------------------------
    // Solar Radiation Pressure Implementation
    //--------------------------------------------------------
    
    SolarRadiationPressure::SolarRadiationPressure(double reflectivityCoefficient, 
                                                 double crossSectionalArea,
                                                 double satelliteMass)
        : reflectivityCoefficient_(reflectivityCoefficient),
          crossSectionalArea_(crossSectionalArea),
          satelliteMass_(satelliteMass) {
        
        if (satelliteMass_ <= 0.0) {
            throw std::invalid_argument("Satellite mass must be positive");
        }
        
        if (crossSectionalArea_ <= 0.0) {
            throw std::invalid_argument("Cross-sectional area must be positive");
        }
        
        if (reflectivityCoefficient_ < 0.0) {
            throw std::invalid_argument("Reflectivity coefficient cannot be negative");
        }
    }
    
    Vector3 SolarRadiationPressure::calculateAcceleration(const CartesianState& state, double time) const {
        if (!enabled_) {
            return Vector3::Zero();
        }
        
        // Get position
        Vector3 position = state.getPosition();
        
        // Calculate Sun position
        Vector3 sunPosition = calculateSunPosition(time);
        
        // Vector from satellite to Sun
        Vector3 satToSun = sunPosition - position;
        double satToSunDist = satToSun.norm();
        
        if (satToSunDist < 1e-10) {
            return Vector3::Zero(); // Avoid division by zero (extremely unlikely)
        }
        
        // Normalize vector to get direction from satellite to Sun
        Vector3 sunDirection = satToSun / satToSunDist;
        
        // Check if satellite is in Earth's shadow
        if (isInShadow(position, sunPosition)) {
            return Vector3::Zero(); // No radiation pressure in shadow
        }
        
        // Solar radiation flux at Earth's distance (W/m²)
        constexpr double solarFlux = 1366.0;
        
        // Speed of light (m/s)
        constexpr double speedOfLight = 299792458.0;
        
        // Calculate solar radiation pressure
        // P = F/c * (1 + reflectivity)
        double pressure = solarFlux / speedOfLight * (1.0 + reflectivityCoefficient_);
        
        // Calculate acceleration: a = F/m = P*A/m
        double accelerationMagnitude = pressure * crossSectionalArea_ / satelliteMass_;
        
        // Return acceleration vector (in the direction from Sun to satellite)
        return -accelerationMagnitude * sunDirection; // Negative because pressure pushes away from Sun
    }
    
    Vector3 SolarRadiationPressure::calculateSunPosition(double time) const {
        // Simplified Sun position model based on time
        // This is a very basic model that assumes the Sun moves in a circular orbit around Earth (which is not true)
        // For a real simulation, an ephemeris model would be used
        
        // Earth-Sun distance (m)
        constexpr double earthSunDistance = 149597870700.0; // 1 AU
        
        // Approximate sidereal year (s)
        constexpr double siderealYear = 365.25636 * 24.0 * 3600.0;
        
        // Angular position of the Sun (simplified circular motion)
        double theta = 2.0 * M_PI * std::fmod(time, siderealYear) / siderealYear;
        
        // Sun position in Earth-centered inertial frame
        double x = earthSunDistance * std::cos(theta);
        double y = earthSunDistance * std::sin(theta);
        double z = 0.0; // Simplified: Sun in the ecliptic plane
        
        return Vector3(x, y, z);
    }
    
    bool SolarRadiationPressure::isInShadow(const Vector3& satellitePosition, const Vector3& sunPosition) const {
        // Vector from Earth to satellite
        Vector3 earthToSat = satellitePosition; // Since Earth is at origin
        double earthToSatDist = earthToSat.norm();
        
        // Vector from Earth to Sun
        Vector3 earthToSun = sunPosition;
        double earthToSunDist = earthToSun.norm();
        
        // Calculate the angle between Earth-to-Satellite and Earth-to-Sun vectors
        double cosAngle = earthToSat.dot(earthToSun) / (earthToSatDist * earthToSunDist);
        
        // If the angle is greater than 90 degrees, satellite is on the opposite side of Earth from Sun
        if (cosAngle < 0) {
            // Check if satellite is within Earth's shadow cone
            double sinUmbra = EARTH_RADIUS / earthToSunDist;
            double sinSat = std::sqrt(1 - cosAngle * cosAngle);
            
            // If sinSat < sinUmbra and satellite is on opposite side, it's in shadow
            return (sinSat < sinUmbra);
        }
        
        return false;
    }
    
    //--------------------------------------------------------
    // Third-Body Perturbation Implementation
    //--------------------------------------------------------
    
    ThirdBodyPerturbation::ThirdBodyPerturbation(Body body)
        : bodyType_(body) {
        
        // Set mass based on body type
        switch (bodyType_) {
            case Body::MOON:
                bodyMass_ = 7.342e22; // Moon mass in kg
                break;
            case Body::SUN:
                bodyMass_ = 1.989e30; // Sun mass in kg
                break;
            default:
                throw std::invalid_argument("Unknown third body type");
        }
    }
    
    std::string ThirdBodyPerturbation::getName() const {
        switch (bodyType_) {
            case Body::MOON:
                return "Moon Perturbation";
            case Body::SUN:
                return "Sun Perturbation";
            default:
                return "Unknown Third-Body Perturbation";
        }
    }
    
    Vector3 ThirdBodyPerturbation::calculateAcceleration(const CartesianState& state, double time) const {
        if (!enabled_) {
            return Vector3::Zero();
        }
        
        // Get satellite position
        Vector3 satellitePosition = state.getPosition();
        
        // Get third body position
        Vector3 bodyPosition = calculateBodyPosition(time);
        
        // Vector from satellite to third body
        Vector3 satToBody = bodyPosition - satellitePosition;
        double satToBodyDist = satToBody.norm();
        
        if (satToBodyDist < 1e-10) {
            return Vector3::Zero(); // Avoid division by zero (extremely unlikely)
        }
        
        // Vector from Earth to third body
        Vector3 earthToBody = bodyPosition;
        double earthToBodyDist = earthToBody.norm();
        
        if (earthToBodyDist < 1e-10) {
            return Vector3::Zero(); // Avoid division by zero (extremely unlikely)
        }
        
        // Gravitational constant
        constexpr double G = 6.67430e-11; // m³/(kg·s²)
        
        // Calculate third-body perturbation acceleration
        // a = GMbody * [ (r_body - r_sat) / |r_body - r_sat|³ - r_body / |r_body|³ ]
        double bodyFactor = G * bodyMass_;
        
        Vector3 term1 = satToBody / std::pow(satToBodyDist, 3);
        Vector3 term2 = earthToBody / std::pow(earthToBodyDist, 3);
        
        return bodyFactor * (term1 - term2);
    }
    
    Vector3 ThirdBodyPerturbation::calculateBodyPosition(double time) const {
        switch (bodyType_) {
            case Body::MOON: {
                // Simplified Moon position model
                // Assumes circular orbit with 27.32 day period at 384,400 km
                
                double moonOrbitalRadius = 384400000.0; // m
                double moonOrbitalPeriod = 27.32 * 24.0 * 3600.0; // s
                
                double theta = 2.0 * M_PI * std::fmod(time, moonOrbitalPeriod) / moonOrbitalPeriod;
                
                // Add 5.145° inclination
                double inclination = 5.145 * M_PI / 180.0; // radians
                
                // Moon position in Earth-centered inertial frame
                double x = moonOrbitalRadius * std::cos(theta);
                double y = moonOrbitalRadius * std::sin(theta) * std::cos(inclination);
                double z = moonOrbitalRadius * std::sin(theta) * std::sin(inclination);
                
                return Vector3(x, y, z);
            }
            
            case Body::SUN: {
                // Re-use the Sun position calculation from the SRP model
                return SolarRadiationPressure().calculateSunPosition(time);
            }
                
            default:
                return Vector3::Zero();
        }
    }

} // namespace SatelliteSimulator