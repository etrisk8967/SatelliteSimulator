// physics_engine/src/orbit_propagator.cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdexcept>
#include <spdlog/spdlog.h>
#include "../include/orbit_propagator.h"

namespace SatelliteSimulator {

    // Constructor with default values
    OrbitPropagator::OrbitPropagator()
        : currentState_(),
          stepSize_(10.0),            // Default 10 seconds
          absoluteTolerance_(1e-10),  // Default absolute tolerance
          relativeTolerance_(1e-8)    // Default relative tolerance
    {
    }

    // Set the initial state
    void OrbitPropagator::setInitialState(const CartesianState& initialState) {
        currentState_ = initialState;
    }

    // Get the current state
    CartesianState OrbitPropagator::getCurrentState() const {
        return currentState_;
    }

    // Set propagation parameters
    void OrbitPropagator::setStepSize(double stepSize) {
        if (stepSize <= 0.0) {
            throw std::invalid_argument("Step size must be positive");
        }
        stepSize_ = stepSize;
    }

    void OrbitPropagator::setAbsoluteTolerance(double absTol) {
        if (absTol <= 0.0) {
            throw std::invalid_argument("Absolute tolerance must be positive");
        }
        absoluteTolerance_ = absTol;
    }

    void OrbitPropagator::setRelativeTolerance(double relTol) {
        if (relTol <= 0.0) {
            throw std::invalid_argument("Relative tolerance must be positive");
        }
        relativeTolerance_ = relTol;
    }

    // Perturbation management
    void OrbitPropagator::addPerturbation(std::shared_ptr<PerturbationForce> perturbation) {
        perturbations_.push_back(perturbation);
    }

    void OrbitPropagator::clearPerturbations() {
        perturbations_.clear();
    }

    std::vector<std::shared_ptr<PerturbationForce>> OrbitPropagator::getPerturbations() const {
        return perturbations_;
    }

    // Maneuver management
    void OrbitPropagator::addManeuver(std::shared_ptr<OrbitManeuver> maneuver) {
        maneuvers_.push_back(maneuver);
    }

    void OrbitPropagator::clearManeuvers() {
        maneuvers_.clear();
    }

    std::vector<std::shared_ptr<OrbitManeuver>> OrbitPropagator::getManeuvers() const {
        return maneuvers_;
    }

    // Execute an impulsive maneuver immediately
    void OrbitPropagator::executeImpulsiveManeuver(const ImpulsiveManeuver& maneuver) {
        currentState_ = maneuver.apply(currentState_);
    }

    // Calculate acceleration due to gravity (two-body model)
    Vector3 OrbitPropagator::calculateGravitationalAcceleration(const Vector3& position, double mu) {
        double r = position.norm();
        
        if (r < 1e-10) {
            throw std::runtime_error("Position too close to origin, gravitational acceleration undefined");
        }
        
        // a = -μ/r³ * r (Newton's law of gravitation)
        return -mu / (r * r * r) * position;
    }

    // Calculate total acceleration (gravity + perturbations + active maneuvers)
    Vector3 OrbitPropagator::calculateTotalAcceleration(const CartesianState& state, double time) const {
        // Start with the basic gravitational acceleration
        Vector3 totalAcceleration = calculateGravitationalAcceleration(state.getPosition());
        
        // Add each enabled perturbation
        for (const auto& perturbation : perturbations_) {
            if (perturbation->isEnabled()) {
                totalAcceleration += perturbation->calculateAcceleration(state, time);
            }
        }
        
        // Add acceleration from active maneuvers (finite burns)
        for (const auto& maneuver : maneuvers_) {
            if (maneuver->isEnabled() && maneuver->isActive(time)) {
                totalAcceleration += maneuver->getAcceleration(state, time);
            }
        }
        
        return totalAcceleration;
    }

    // Convert between CartesianState and state vector for integrator
    // State vector format: [px, py, pz, vx, vy, vz]
    StateVector OrbitPropagator::cartesianToStateVector(const CartesianState& state) {
        Vector3 position = state.getPosition();
        Vector3 velocity = state.getVelocity();
        
        StateVector stateVector(6);
        stateVector[0] = position.x();
        stateVector[1] = position.y();
        stateVector[2] = position.z();
        stateVector[3] = velocity.x();
        stateVector[4] = velocity.y();
        stateVector[5] = velocity.z();
        
        return stateVector;
    }

    CartesianState OrbitPropagator::stateVectorToCartesian(const StateVector& stateVector) {
        if (stateVector.size() != 6) {
            throw std::invalid_argument("State vector must have size 6");
        }
        
        Vector3 position(stateVector[0], stateVector[1], stateVector[2]);
        Vector3 velocity(stateVector[3], stateVector[4], stateVector[5]);
        
        return CartesianState(position, velocity);
    }

    // Define the system of differential equations
    // x' = f(x, t) where x = [px, py, pz, vx, vy, vz]
    // For two-body problem with perturbations and maneuvers:
    // p' = v
    // v' = a_gravity + a_perturbations + a_maneuvers
    void OrbitPropagator::systemFunction(const StateVector& x, StateVector& dxdt, double t) const {
        // Convert state vector to CartesianState
        CartesianState state = stateVectorToCartesian(x);
        
        // Calculate total acceleration including perturbations and maneuvers
        Vector3 acceleration = calculateTotalAcceleration(state, t);
        
        // Set derivatives
        dxdt[0] = x[3];  // dx/dt = vx
        dxdt[1] = x[4];  // dy/dt = vy
        dxdt[2] = x[5];  // dz/dt = vz
        dxdt[3] = acceleration.x();  // dvx/dt = ax
        dxdt[4] = acceleration.y();  // dvy/dt = ay
        dxdt[5] = acceleration.z();  // dvz/dt = az
    }

    // Check for and execute any impulsive maneuvers at the specified time
    void OrbitPropagator::checkAndExecuteImpulsiveManeuvers(double time) {
        for (const auto& maneuver : maneuvers_) {
            // Check if this is an impulsive maneuver and if it's time to execute it
            auto impulsiveManeuver = std::dynamic_pointer_cast<ImpulsiveManeuver>(maneuver);
            if (impulsiveManeuver && impulsiveManeuver->isEnabled() && 
                std::abs(time - impulsiveManeuver->getStartTime()) < 1e-10) {
                
                // Execute the impulsive maneuver
                currentState_ = impulsiveManeuver->apply(currentState_);
                
                // Log the maneuver execution
                spdlog::info("Executed impulsive maneuver at time {:.2f}s: {:.2f} m/s {}",
                            time, impulsiveManeuver->getDeltaV(), impulsiveManeuver->getName());
            }
        }
    }

    // Propagate the orbital state to a specific time
    CartesianState OrbitPropagator::propagateToTime(double finalTime, const StateObserver& observer) {
        return propagate(0.0, finalTime, observer);
    }

    // Propagate the orbital state by a time interval
    CartesianState OrbitPropagator::propagateByTime(double deltaTime, const StateObserver& observer) {
        return propagate(0.0, deltaTime, observer);
    }

    // Core propagation implementation
    CartesianState OrbitPropagator::propagate(double startTime, double endTime, const StateObserver& observer) {
        if (endTime < startTime) {
            throw std::invalid_argument("End time must be greater than or equal to start time");
        }
        
        if (std::abs(endTime - startTime) < 1e-10) {
            // No propagation needed
            return currentState_;
        }
        
        // Check if we have any impulsive maneuvers that need to be executed during this propagation
        std::vector<double> maneuverTimes;
        for (const auto& maneuver : maneuvers_) {
            auto impulsiveManeuver = std::dynamic_pointer_cast<ImpulsiveManeuver>(maneuver);
            if (impulsiveManeuver && impulsiveManeuver->isEnabled()) {
                double mTime = impulsiveManeuver->getStartTime();
                if (mTime >= startTime && mTime <= endTime) {
                    maneuverTimes.push_back(mTime);
                }
            }
        }
        
        // If we have impulsive maneuvers, we need to propagate in segments
        if (!maneuverTimes.empty()) {
            // Sort maneuver times
            std::sort(maneuverTimes.begin(), maneuverTimes.end());
            
            // Start with current time
            double currentTime = startTime;
            
            // Propagate to each maneuver time, execute the maneuver, and continue
            for (double maneuverTime : maneuverTimes) {
                // Propagate to the maneuver time
                if (maneuverTime > currentTime) {
                    // Create an observer that adjusts the time
                    auto timeAdjustedObserver = observer ? 
                        [&observer, currentTime](const CartesianState& state, double t) {
                            observer(state, t + currentTime);
                        } : StateObserver(nullptr);
                    
                    // Propagate to the maneuver time
                    propagate(0.0, maneuverTime - currentTime, timeAdjustedObserver);
                }
                
                // Execute any impulsive maneuvers at this time
                checkAndExecuteImpulsiveManeuvers(maneuverTime);
                
                // If observer is provided, call it at the maneuver time
                if (observer) {
                    observer(currentState_, maneuverTime);
                }
                
                // Update current time
                currentTime = maneuverTime;
            }
            
            // Propagate to the end time
            if (endTime > currentTime) {
                // Create an observer that adjusts the time
                auto timeAdjustedObserver = observer ? 
                    [&observer, currentTime](const CartesianState& state, double t) {
                        observer(state, t + currentTime);
                    } : StateObserver(nullptr);
                
                // Propagate to the end time
                propagate(0.0, endTime - currentTime, timeAdjustedObserver);
            }
            
            return currentState_;
        }
        
        // If we don't have any impulsive maneuvers, propagate normally
        
        // Convert current state to state vector
        StateVector state = cartesianToStateVector(currentState_);
        
        // Create stepper - use controlled Dormand-Prince 5 (DOPRI5) with error control
        using namespace boost::numeric::odeint;
        
        auto stepper = make_controlled(absoluteTolerance_, relativeTolerance_, 
                                      runge_kutta_dopri5<StateVector>());
        
        // Create system function wrapper
        SystemFunctionWrapper system(*this);
        
        // Lambda to observe state during integration
        auto observerFunction = [this, &observer](const StateVector& x, double t) {
            if (observer) {
                CartesianState state = stateVectorToCartesian(x);
                observer(state, t);
            }
        };
        
        // Propagate the orbit
        if (observer) {
            integrate_adaptive(stepper, system, state, startTime, endTime, stepSize_, observerFunction);
        } else {
            integrate_adaptive(stepper, system, state, startTime, endTime, stepSize_);
        }
        
        // Convert back to CartesianState
        currentState_ = stateVectorToCartesian(state);
        
        return currentState_;
    }

} // namespace SatelliteSimulator