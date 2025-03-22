// physics_engine/include/orbit_maneuver.h
#pragma once

#include <string>
#include <memory>
#include <functional>
#include "vector_math.h"
#include "orbital_state.h"

namespace SatelliteSimulator {

    // Forward declaration
    class OrbitPropagator;

    // Base class for all orbit maneuvers
    class OrbitManeuver {
    public:
        // Maneuver direction options
        enum class Direction {
            PROGRADE,       // Along velocity vector
            RETROGRADE,     // Opposite to velocity vector
            NORMAL,         // Perpendicular to orbital plane (northward)
            ANTI_NORMAL,    // Perpendicular to orbital plane (southward)
            RADIAL_IN,      // Toward central body
            RADIAL_OUT,     // Away from central body
            CUSTOM          // Custom direction specified by a vector
        };
        
        virtual ~OrbitManeuver() = default;
        
        // Pure virtual methods to be implemented by derived classes
        virtual Vector3 getAcceleration(const CartesianState& state, double time) const = 0;
        virtual double getFuelConsumed(double time) const = 0;
        virtual bool isActive(double time) const = 0;
        virtual std::string getName() const = 0;
        
        // Getter/Setter for execution time
        double getStartTime() const { return startTime_; }
        void setStartTime(double startTime) { startTime_ = startTime; }
        
        // Enable/disable maneuver
        bool isEnabled() const { return enabled_; }
        void setEnabled(bool enabled) { enabled_ = enabled; }
        
    protected:
        // Constructor is protected as this is an abstract base class
        OrbitManeuver(double startTime = 0.0)
            : startTime_(startTime), enabled_(true) {}
        
        // Common utility to calculate direction vector based on maneuver direction
        Vector3 calculateDirectionVector(const CartesianState& state, Direction direction, 
                                         const Vector3& customDirection = Vector3::Zero()) const;
        
        double startTime_;  // Time at which the maneuver begins (seconds)
        bool enabled_;      // Whether the maneuver is enabled
    };

    // Impulsive maneuver (instantaneous change in velocity)
    class ImpulsiveManeuver : public OrbitManeuver {
    public:
        ImpulsiveManeuver(double deltaV, Direction direction, 
                         double startTime = 0.0, 
                         const Vector3& customDirection = Vector3::Zero());
        
        // Implementation of base class virtual methods
        Vector3 getAcceleration(const CartesianState& state, double time) const override;
        double getFuelConsumed(double time) const override;
        bool isActive(double time) const override;
        std::string getName() const override { return "Impulsive Maneuver"; }
        
        // Getters/Setters
        double getDeltaV() const { return deltaV_; }
        void setDeltaV(double deltaV) { deltaV_ = deltaV; }
        
        Direction getDirection() const { return direction_; }
        void setDirection(Direction direction) { direction_ = direction; }
        
        const Vector3& getCustomDirection() const { return customDirection_; }
        void setCustomDirection(const Vector3& direction) { customDirection_ = direction; }
        
        // Apply the maneuver directly to a state
        CartesianState apply(const CartesianState& state) const;
        
    private:
        double deltaV_;             // Magnitude of velocity change (m/s)
        Direction direction_;       // Direction of the maneuver
        Vector3 customDirection_;   // Custom direction vector (if Direction::CUSTOM)
    };

    // Finite burn maneuver (continuous thrust over time)
    class FiniteBurnManeuver : public OrbitManeuver {
    public:
        FiniteBurnManeuver(double thrust, double specificImpulse, double satelliteMass,
                          double duration, Direction direction,
                          double startTime = 0.0,
                          const Vector3& customDirection = Vector3::Zero());
        
        // Implementation of base class virtual methods
        Vector3 getAcceleration(const CartesianState& state, double time) const override;
        double getFuelConsumed(double time) const override;
        bool isActive(double time) const override;
        std::string getName() const override { return "Finite Burn Maneuver"; }
        
        // Getters/Setters
        double getThrust() const { return thrust_; }
        void setThrust(double thrust) { thrust_ = thrust; }
        
        double getSpecificImpulse() const { return specificImpulse_; }
        void setSpecificImpulse(double isp) { specificImpulse_ = isp; }
        
        double getSatelliteMass() const { return satelliteMass_; }
        void setSatelliteMass(double mass) { satelliteMass_ = mass; }
        
        double getDuration() const { return duration_; }
        void setDuration(double duration) { duration_ = duration; }
        
        Direction getDirection() const { return direction_; }
        void setDirection(Direction direction) { direction_ = direction; }
        
        const Vector3& getCustomDirection() const { return customDirection_; }
        void setCustomDirection(const Vector3& direction) { customDirection_ = direction; }
        
        // Calculate the total delta-V that will be provided by this burn
        double calculateTotalDeltaV() const;
        
        // Calculate the current mass at a given time (accounting for propellant usage)
        double calculateCurrentMass(double time) const;
        
    private:
        double thrust_;             // Thrust force (N)
        double specificImpulse_;    // Specific impulse (s)
        double satelliteMass_;      // Initial satellite mass (kg)
        double duration_;           // Duration of the burn (s)
        Direction direction_;       // Direction of the maneuver
        Vector3 customDirection_;   // Custom direction vector (if Direction::CUSTOM)
        
        // Constants
        static constexpr double GRAVITY = 9.80665;  // Standard gravity (m/s²)
    };

    // Utility functions for common maneuver calculations
    
    // Calculate delta-V required for a Hohmann transfer between circular orbits
    double calculateHohmannTransferDeltaV(double initialRadius, double finalRadius, double mu = EARTH_MU);
    
    // Calculate delta-V required for a plane change
    double calculatePlaneChangeDeltaV(double velocity, double planeAngle);
    
    // Calculate delta-V required for a combined maneuver (plane change during Hohmann transfer)
    double calculateCombinedManeuverDeltaV(double initialRadius, double finalRadius, 
                                          double planeAngle, double mu = EARTH_MU);
    
    // Create an impulsive maneuver for Hohmann transfer (first burn)
    std::shared_ptr<ImpulsiveManeuver> createHohmannTransferInitialBurn(
        double initialRadius, double finalRadius, double startTime = 0.0, double mu = EARTH_MU);
    
    // Create an impulsive maneuver for Hohmann transfer (second burn)
    std::shared_ptr<ImpulsiveManeuver> createHohmannTransferFinalBurn(
        double initialRadius, double finalRadius, double startTime = 0.0, double mu = EARTH_MU);
    
    // Create an impulsive maneuver for plane change
    std::shared_ptr<ImpulsiveManeuver> createPlaneChangeManeuver(
        double velocity, double planeAngle, double startTime = 0.0);
    
    // Calculate time of flight for a Hohmann transfer
    double calculateHohmannTimeOfFlight(double initialRadius, double finalRadius, double mu = EARTH_MU);

} // namespace SatelliteSimulator