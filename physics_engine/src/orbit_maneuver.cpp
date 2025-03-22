// physics_engine/src/orbit_maneuver.cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdexcept>
#include <spdlog/spdlog.h>
#include "../include/orbit_maneuver.h"

namespace SatelliteSimulator {

    //--------------------------------------------------------
    // OrbitManeuver Base Class
    //--------------------------------------------------------
    
    Vector3 OrbitManeuver::calculateDirectionVector(const CartesianState& state, Direction direction, 
                                                   const Vector3& customDirection) const {
        // Get position and velocity vectors
        Vector3 position = state.getPosition();
        Vector3 velocity = state.getVelocity();
        
        // Ensure we have non-zero vectors
        double positionNorm = position.norm();
        double velocityNorm = velocity.norm();
        
        if (positionNorm < 1e-10 || velocityNorm < 1e-10) {
            throw std::runtime_error("Position or velocity too close to zero for direction calculation");
        }
        
        // Calculate the orbital frame unit vectors
        Vector3 vHat = velocity / velocityNorm;            // Velocity direction (prograde)
        Vector3 hVec = position.cross(velocity);           // Angular momentum vector
        double hNorm = hVec.norm();
        
        if (hNorm < 1e-10) {
            throw std::runtime_error("Angular momentum too close to zero for direction calculation");
        }
        
        Vector3 hHat = hVec / hNorm;                      // Angular momentum direction (normal)
        Vector3 rHat = position / positionNorm;           // Radial direction (outward)
        
        // Return unit vector in the requested direction
        switch (direction) {
            case Direction::PROGRADE:
                return vHat;
            
            case Direction::RETROGRADE:
                return -vHat;
            
            case Direction::NORMAL:
                return hHat;
            
            case Direction::ANTI_NORMAL:
                return -hHat;
            
            case Direction::RADIAL_OUT:
                return rHat;
            
            case Direction::RADIAL_IN:
                return -rHat;
            
            case Direction::CUSTOM:
                // Normalize custom direction if provided
                if (customDirection.norm() < 1e-10) {
                    throw std::runtime_error("Custom direction vector too close to zero");
                }
                return customDirection.normalized();
            
            default:
                throw std::runtime_error("Unknown maneuver direction");
        }
    }

    //--------------------------------------------------------
    // ImpulsiveManeuver Implementation
    //--------------------------------------------------------
    
    ImpulsiveManeuver::ImpulsiveManeuver(double deltaV, Direction direction, 
                                       double startTime, 
                                       const Vector3& customDirection)
        : OrbitManeuver(startTime),
          deltaV_(deltaV),
          direction_(direction),
          customDirection_(customDirection) {
        
        if (deltaV_ < 0.0) {
            throw std::invalid_argument("Delta-V cannot be negative");
        }
        
        if (direction_ == Direction::CUSTOM && customDirection_.norm() < 1e-10) {
            throw std::invalid_argument("Custom direction vector too close to zero");
        }
    }
    
    Vector3 ImpulsiveManeuver::getAcceleration(const CartesianState& state, double time) const {
        // Impulsive maneuvers are instantaneous, so they don't provide continuous acceleration
        // This method should not be called for impulsive maneuvers
        return Vector3::Zero();
    }
    
    double ImpulsiveManeuver::getFuelConsumed(double time) const {
        // Fuel consumption is instantaneous at the moment of the burn
        // For simplicity, we don't calculate fuel consumption for impulsive maneuvers
        return 0.0;
    }
    
    bool ImpulsiveManeuver::isActive(double time) const {
        // Impulsive maneuvers are only active at exactly the start time
        // In practice, they're applied directly to the state and not via continuous acceleration
        return false;
    }
    
    CartesianState ImpulsiveManeuver::apply(const CartesianState& state) const {
        if (!enabled_) {
            return state;  // No change if maneuver is disabled
        }
        
        // Calculate the direction vector for the maneuver
        Vector3 directionVector = calculateDirectionVector(state, direction_, customDirection_);
        
        // Apply the velocity change
        Vector3 newVelocity = state.getVelocity() + directionVector * deltaV_;
        
        // Return new state with the same position but updated velocity
        return CartesianState(state.getPosition(), newVelocity);
    }

    //--------------------------------------------------------
    // FiniteBurnManeuver Implementation
    //--------------------------------------------------------
    
    FiniteBurnManeuver::FiniteBurnManeuver(double thrust, double specificImpulse, double satelliteMass,
                                         double duration, Direction direction,
                                         double startTime,
                                         const Vector3& customDirection)
        : OrbitManeuver(startTime),
          thrust_(thrust),
          specificImpulse_(specificImpulse),
          satelliteMass_(satelliteMass),
          duration_(duration),
          direction_(direction),
          customDirection_(customDirection) {
        
        if (thrust_ < 0.0) {
            throw std::invalid_argument("Thrust cannot be negative");
        }
        
        if (specificImpulse_ <= 0.0) {
            throw std::invalid_argument("Specific impulse must be positive");
        }
        
        if (satelliteMass_ <= 0.0) {
            throw std::invalid_argument("Satellite mass must be positive");
        }
        
        if (duration_ <= 0.0) {
            throw std::invalid_argument("Duration must be positive");
        }
        
        if (direction_ == Direction::CUSTOM && customDirection_.norm() < 1e-10) {
            throw std::invalid_argument("Custom direction vector too close to zero");
        }
    }
    
    Vector3 FiniteBurnManeuver::getAcceleration(const CartesianState& state, double time) const {
        if (!enabled_ || !isActive(time)) {
            return Vector3::Zero();  // No acceleration if not active
        }
        
        // Calculate the direction vector for the maneuver
        Vector3 directionVector = calculateDirectionVector(state, direction_, customDirection_);
        
        // Calculate the current mass
        double currentMass = calculateCurrentMass(time);
        
        // Calculate the acceleration magnitude (F = m*a => a = F/m)
        double accelerationMagnitude = thrust_ / currentMass;
        
        // Return the acceleration vector
        return directionVector * accelerationMagnitude;
    }
    
    double FiniteBurnManeuver::getFuelConsumed(double time) const {
        if (!enabled_) {
            return 0.0;  // No fuel consumed if disabled
        }
        
        if (time <= startTime_) {
            return 0.0;  // No fuel consumed before the maneuver starts
        }
        
        double endTime = startTime_ + duration_;
        double burnTime = std::min(time, endTime) - startTime_;
        
        // Calculate mass flow rate (kg/s)
        double massFlowRate = thrust_ / (specificImpulse_ * GRAVITY);
        
        // Return total fuel consumed (kg)
        return massFlowRate * burnTime;
    }
    
    bool FiniteBurnManeuver::isActive(double time) const {
        if (!enabled_) {
            return false;
        }
        
        return (time >= startTime_ && time < startTime_ + duration_);
    }
    
    double FiniteBurnManeuver::calculateTotalDeltaV() const {
        // Calculate the delta-V using the rocket equation: ΔV = Isp * g * ln(m₀ / mf)
        double massFlowRate = thrust_ / (specificImpulse_ * GRAVITY);
        double finalMass = satelliteMass_ - massFlowRate * duration_;
        
        if (finalMass <= 0.0) {
            throw std::runtime_error("Burn depletes all fuel and more, recalculate duration");
        }
        
        return specificImpulse_ * GRAVITY * std::log(satelliteMass_ / finalMass);
    }
    
    double FiniteBurnManeuver::calculateCurrentMass(double time) const {
        if (!enabled_ || time <= startTime_) {
            return satelliteMass_;  // Initial mass before the burn
        }
        
        double endTime = startTime_ + duration_;
        if (time >= endTime) {
            // Calculate final mass after the complete burn
            double massFlowRate = thrust_ / (specificImpulse_ * GRAVITY);
            return satelliteMass_ - massFlowRate * duration_;
        }
        
        // Calculate current mass during the burn
        double burnTime = time - startTime_;
        double massFlowRate = thrust_ / (specificImpulse_ * GRAVITY);
        return satelliteMass_ - massFlowRate * burnTime;
    }

    //--------------------------------------------------------
    // Utility Functions
    //--------------------------------------------------------
    
    double calculateHohmannTransferDeltaV(double initialRadius, double finalRadius, double mu) {
        if (initialRadius <= 0.0 || finalRadius <= 0.0) {
            throw std::invalid_argument("Orbit radii must be positive");
        }
        
        // Calculate the semi-major axis of the transfer orbit
        double transferSMA = (initialRadius + finalRadius) / 2.0;
        
        // Calculate velocities in the initial circular orbit
        double initialVelocity = std::sqrt(mu / initialRadius);
        
        // Calculate velocities at periapsis of the transfer orbit
        double transferVelocityAtPeriapsis = std::sqrt(mu * (2.0 / initialRadius - 1.0 / transferSMA));
        
        // Calculate velocities at apoapsis of the transfer orbit
        double transferVelocityAtApoapsis = std::sqrt(mu * (2.0 / finalRadius - 1.0 / transferSMA));
        
        // Calculate velocities in the final circular orbit
        double finalVelocity = std::sqrt(mu / finalRadius);
        
        // Calculate delta-V for first burn (at periapsis of transfer orbit)
        double deltaV1 = std::abs(transferVelocityAtPeriapsis - initialVelocity);
        
        // Calculate delta-V for second burn (at apoapsis of transfer orbit)
        double deltaV2 = std::abs(finalVelocity - transferVelocityAtApoapsis);
        
        // Return total delta-V
        return deltaV1 + deltaV2;
    }
    
    double calculatePlaneChangeDeltaV(double velocity, double planeAngle) {
        if (velocity < 0.0) {
            throw std::invalid_argument("Velocity must be non-negative");
        }
        
        // Delta-V for a plane change = 2 * v * sin(angle/2)
        return 2.0 * velocity * std::sin(planeAngle / 2.0);
    }
    
    double calculateCombinedManeuverDeltaV(double initialRadius, double finalRadius, 
                                          double planeAngle, double mu) {
        if (initialRadius <= 0.0 || finalRadius <= 0.0) {
            throw std::invalid_argument("Orbit radii must be positive");
        }
        
        // Calculate the semi-major axis of the transfer orbit
        double transferSMA = (initialRadius + finalRadius) / 2.0;
        
        // Calculate velocities in the initial circular orbit
        double initialVelocity = std::sqrt(mu / initialRadius);
        
        // Calculate velocities at periapsis of the transfer orbit
        double transferVelocityAtPeriapsis = std::sqrt(mu * (2.0 / initialRadius - 1.0 / transferSMA));
        
        // Calculate velocities at apoapsis of the transfer orbit
        double transferVelocityAtApoapsis = std::sqrt(mu * (2.0 / finalRadius - 1.0 / transferSMA));
        
        // Calculate velocities in the final circular orbit
        double finalVelocity = std::sqrt(mu / finalRadius);
        
        // Calculate delta-V for first burn (at periapsis of transfer orbit)
        double deltaV1 = std::abs(transferVelocityAtPeriapsis - initialVelocity);
        
        // Calculate delta-V for combined second burn and plane change
        // For a combined maneuver, we use: ΔV = √(v₁² + v₂² - 2*v₁*v₂*cos(angle))
        double deltaV2 = std::sqrt(transferVelocityAtApoapsis * transferVelocityAtApoapsis + 
                                  finalVelocity * finalVelocity - 
                                  2.0 * transferVelocityAtApoapsis * finalVelocity * std::cos(planeAngle));
        
        // Return total delta-V
        return deltaV1 + deltaV2;
    }
    
    std::shared_ptr<ImpulsiveManeuver> createHohmannTransferInitialBurn(
        double initialRadius, double finalRadius, double startTime, double mu) {
        
        if (initialRadius <= 0.0 || finalRadius <= 0.0) {
            throw std::invalid_argument("Orbit radii must be positive");
        }
        
        // Calculate the semi-major axis of the transfer orbit
        double transferSMA = (initialRadius + finalRadius) / 2.0;
        
        // Calculate velocities in the initial circular orbit
        double initialVelocity = std::sqrt(mu / initialRadius);
        
        // Calculate velocities at periapsis of the transfer orbit
        double transferVelocityAtPeriapsis = std::sqrt(mu * (2.0 / initialRadius - 1.0 / transferSMA));
        
        // Calculate delta-V for first burn
        double deltaV = transferVelocityAtPeriapsis - initialVelocity;
        
        // Create the maneuver with the appropriate direction
        OrbitManeuver::Direction direction = (deltaV >= 0.0) ? 
            OrbitManeuver::Direction::PROGRADE : OrbitManeuver::Direction::RETROGRADE;
        
        return std::make_shared<ImpulsiveManeuver>(std::abs(deltaV), direction, startTime);
    }
    
    std::shared_ptr<ImpulsiveManeuver> createHohmannTransferFinalBurn(
        double initialRadius, double finalRadius, double startTime, double mu) {
        
        if (initialRadius <= 0.0 || finalRadius <= 0.0) {
            throw std::invalid_argument("Orbit radii must be positive");
        }
        
        // Calculate the semi-major axis of the transfer orbit
        double transferSMA = (initialRadius + finalRadius) / 2.0;
        
        // Calculate velocities at apoapsis of the transfer orbit
        double transferVelocityAtApoapsis = std::sqrt(mu * (2.0 / finalRadius - 1.0 / transferSMA));
        
        // Calculate velocities in the final circular orbit
        double finalVelocity = std::sqrt(mu / finalRadius);
        
        // Calculate delta-V for second burn
        double deltaV = finalVelocity - transferVelocityAtApoapsis;
        
        // Create the maneuver with the appropriate direction
        OrbitManeuver::Direction direction = (deltaV >= 0.0) ? 
            OrbitManeuver::Direction::PROGRADE : OrbitManeuver::Direction::RETROGRADE;
        
        return std::make_shared<ImpulsiveManeuver>(std::abs(deltaV), direction, startTime);
    }
    
    std::shared_ptr<ImpulsiveManeuver> createPlaneChangeManeuver(
        double velocity, double planeAngle, double startTime) {
        
        if (velocity < 0.0) {
            throw std::invalid_argument("Velocity must be non-negative");
        }
        
        // Calculate delta-V for the plane change
        double deltaV = 2.0 * velocity * std::sin(planeAngle / 2.0);
        
        // Create the maneuver (always normal to the orbital plane)
        return std::make_shared<ImpulsiveManeuver>(deltaV, OrbitManeuver::Direction::NORMAL, startTime);
    }
    
    double calculateHohmannTimeOfFlight(double initialRadius, double finalRadius, double mu) {
        if (initialRadius <= 0.0 || finalRadius <= 0.0) {
            throw std::invalid_argument("Orbit radii must be positive");
        }
        
        // Calculate the semi-major axis of the transfer orbit
        double transferSMA = (initialRadius + finalRadius) / 2.0;
        
        // Calculate the time of flight for half an orbit in the transfer ellipse
        // T = π * √(a³/μ)
        return M_PI * std::sqrt(std::pow(transferSMA, 3) / mu);
    }

} // namespace SatelliteSimulator