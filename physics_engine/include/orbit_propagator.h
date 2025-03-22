// physics_engine/include/orbit_propagator.h
#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include "orbital_state.h"
#include "perturbation_force.h"
#include "orbit_maneuver.h"  // Add this include

namespace SatelliteSimulator {

    // Type alias for the state vector
    using StateVector = std::vector<double>;

    // Type alias for observer function that can be called during propagation
    using StateObserver = std::function<void(const CartesianState&, double)>;

    class OrbitPropagator {
    public:
        // Constructor
        OrbitPropagator();

        // Set the initial state
        void setInitialState(const CartesianState& initialState);
        
        // Get the current state
        CartesianState getCurrentState() const;
        
        // Set propagation parameters
        void setStepSize(double stepSize);
        void setAbsoluteTolerance(double absTol);
        void setRelativeTolerance(double relTol);
        
        // Perturbation management
        void addPerturbation(std::shared_ptr<PerturbationForce> perturbation);
        void clearPerturbations();
        std::vector<std::shared_ptr<PerturbationForce>> getPerturbations() const;
        
        // Maneuver management
        void addManeuver(std::shared_ptr<OrbitManeuver> maneuver);
        void clearManeuvers();
        std::vector<std::shared_ptr<OrbitManeuver>> getManeuvers() const;
        
        // Propagation methods
        CartesianState propagateToTime(double finalTime, const StateObserver& observer = nullptr);
        CartesianState propagateByTime(double deltaTime, const StateObserver& observer = nullptr);
        
        // Execute an impulsive maneuver immediately (without propagation)
        void executeImpulsiveManeuver(const ImpulsiveManeuver& maneuver);
        
        // Calculate acceleration due to gravity (two-body model)
        static Vector3 calculateGravitationalAcceleration(const Vector3& position, double mu = EARTH_MU);
        
        // Calculate total acceleration (gravity + perturbations + active maneuvers)
        Vector3 calculateTotalAcceleration(const CartesianState& state, double time) const;
        
    private:
        // Current state
        CartesianState currentState_;
        
        // Integration parameters
        double stepSize_;
        double absoluteTolerance_;
        double relativeTolerance_;
        
        // Perturbation forces
        std::vector<std::shared_ptr<PerturbationForce>> perturbations_;
        
        // Maneuvers
        std::vector<std::shared_ptr<OrbitManeuver>> maneuvers_;
        
        // Convert between CartesianState and state vector for integrator
        static StateVector cartesianToStateVector(const CartesianState& state);
        static CartesianState stateVectorToCartesian(const StateVector& stateVector);
        
        // System function for the integrator
        void systemFunction(const StateVector& x, StateVector& dxdt, double t) const;
        
        // Wrapper for the system function (needed for Boost.Odeint)
        struct SystemFunctionWrapper {
            const OrbitPropagator& propagator;
            
            SystemFunctionWrapper(const OrbitPropagator& prop) : propagator(prop) {}
            
            void operator()(const StateVector& x, StateVector& dxdt, double t) const {
                propagator.systemFunction(x, dxdt, t);
            }
        };
        
        // Check for and execute any impulsive maneuvers at the specified time
        void checkAndExecuteImpulsiveManeuvers(double time);
        
        // Propagate using the integrator
        CartesianState propagate(double startTime, double endTime, const StateObserver& observer);
    };

} // namespace SatelliteSimulator