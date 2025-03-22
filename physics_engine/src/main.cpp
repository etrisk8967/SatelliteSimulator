// physics_engine/src/main.cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include <spdlog/spdlog.h>
#include "../include/vector_math.h"
#include "../include/orbital_state.h"
#include "../include/orbit_propagator.h"
#include "../include/perturbation_force.h"
#include "../include/orbit_maneuver.h"  // New include for orbit maneuvers

using namespace SatelliteSimulator;

void testVectorMath() {
    spdlog::info("Testing Vector Math Utilities");
    
    // Test vector operations
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(4.0, 5.0, 6.0);
    
    spdlog::info("Vector v1: {}", vectorToString(v1));
    spdlog::info("Vector v2: {}", vectorToString(v2));
    
    // Basic operations
    spdlog::info("v1 + v2 = {}", vectorToString(addVectors(v1, v2)));
    spdlog::info("v1 - v2 = {}", vectorToString(subtractVectors(v1, v2)));
    spdlog::info("2 * v1 = {}", vectorToString(scaleVector(v1, 2.0)));
    spdlog::info("v1 • v2 = {}", dotProduct(v1, v2));
    spdlog::info("v1 × v2 = {}", vectorToString(crossProduct(v1, v2)));
    
    // Vector properties
    spdlog::info("Magnitude of v1: {:.6f}", magnitude(v1));
    spdlog::info("Normalized v1: {}", vectorToString(normalize(v1)));
    spdlog::info("Angle between v1 and v2: {:.6f} radians", angleBetweenVectors(v1, v2));
    
    // Coordinate transformations
    spdlog::info("v1 in spherical coordinates: {}", 
                 vectorToString(transformCoordinates(v1, CoordinateSystem::CARTESIAN, CoordinateSystem::SPHERICAL)));
    spdlog::info("v1 in cylindrical coordinates: {}", 
                 vectorToString(transformCoordinates(v1, CoordinateSystem::CARTESIAN, CoordinateSystem::CYLINDRICAL)));
    
    // Test a full coordinate transformation cycle (should get back original vector)
    Vector3 spherical = transformCoordinates(v1, CoordinateSystem::CARTESIAN, CoordinateSystem::SPHERICAL);
    Vector3 backToCartesian = transformCoordinates(spherical, CoordinateSystem::SPHERICAL, CoordinateSystem::CARTESIAN);
    spdlog::info("Original: {} -> Spherical: {} -> Back to Cartesian: {}", 
                 vectorToString(v1), vectorToString(spherical), vectorToString(backToCartesian));
    
    // Rotations
    double angle = M_PI / 4.0;  // 45 degrees
    spdlog::info("Rotating v1 by {} radians around X-axis: {}", 
                 angle, vectorToString(rotateVector(v1, rotationMatrixX(angle))));
    spdlog::info("Rotating v1 by {} radians around Y-axis: {}", 
                 angle, vectorToString(rotateVector(v1, rotationMatrixY(angle))));
    spdlog::info("Rotating v1 by {} radians around Z-axis: {}", 
                 angle, vectorToString(rotateVector(v1, rotationMatrixZ(angle))));
}

void testOrbitalState() {
    spdlog::info("\nTesting Orbital State Representations");
    
    // Test CartesianState
    Vector3 position(7000000.0, 0.0, 0.0);  // 7000 km position along x-axis
    Vector3 velocity(0.0, 7546.0, 0.0);     // Approximate circular orbit velocity along y-axis
    
    CartesianState cartState(position, velocity);
    spdlog::info("Cartesian State:");
    spdlog::info("{}", cartState.toString());
    
    double satelliteMass = 1000.0;  // 1000 kg
    spdlog::info("Kinetic Energy: {:.2f} J", cartState.getKineticEnergy(satelliteMass));
    spdlog::info("Potential Energy: {:.2f} J", cartState.getPotentialEnergy(satelliteMass));
    spdlog::info("Total Energy: {:.2f} J", cartState.getTotalEnergy(satelliteMass));
    spdlog::info("Angular Momentum: {:.2f} kg·m²/s", cartState.getAngularMomentum(satelliteMass));
    
    // Convert to orbital elements
    OrbitalElements elements = cartesianToOrbitalElements(cartState);
    spdlog::info("\nConverted to Orbital Elements:");
    spdlog::info("{}", elements.toString());
    
    // Convert back to Cartesian
    CartesianState convertedCartState = orbitalElementsToCartesian(elements);
    spdlog::info("\nConverted back to Cartesian:");
    spdlog::info("{}", convertedCartState.toString());
    
    // Create orbital elements for an elliptical orbit
    OrbitalElements ellipticalOrbit(8000000.0,    // 8000 km semi-major axis
                                    0.2,          // 0.2 eccentricity
                                    0.5,          // ~28.6 degrees inclination
                                    0.3,          // ~17.2 degrees RAAN
                                    0.4,          // ~22.9 degrees argument of perigee
                                    0.0);         // 0 degrees true anomaly (at periapsis)
    
    spdlog::info("\nElliptical Orbit Elements:");
    spdlog::info("{}", ellipticalOrbit.toString());
    
    // Convert to Cartesian
    CartesianState ellipticalCartState = orbitalElementsToCartesian(ellipticalOrbit);
    spdlog::info("\nElliptical Orbit in Cartesian Coordinates:");
    spdlog::info("{}", ellipticalCartState.toString());
    
    // Test anomaly conversions
    double meanAnomaly = M_PI / 4.0;  // 45 degrees
    double eccentricity = 0.2;
    
    double eccentricAnomaly = meanToEccentricAnomaly(meanAnomaly, eccentricity);
    double trueAnomaly = eccentricToTrueAnomaly(eccentricAnomaly, eccentricity);
    double backToMean = trueToMeanAnomaly(trueAnomaly, eccentricity);
    
    spdlog::info("\nAnomaly Conversions:");
    spdlog::info("Mean Anomaly: {:.6f} rad -> Eccentric Anomaly: {:.6f} rad -> True Anomaly: {:.6f} rad",
                 meanAnomaly, eccentricAnomaly, trueAnomaly);
    spdlog::info("True Anomaly: {:.6f} rad -> Mean Anomaly: {:.6f} rad (original: {:.6f} rad)",
                 trueAnomaly, backToMean, meanAnomaly);
    
    // Test a hyperbolic orbit
    OrbitalElements hyperbolicOrbit(-10000000.0,  // -10000 km semi-major axis (negative for hyperbolic)
                                   1.5,          // 1.5 eccentricity (> 1 for hyperbolic)
                                   0.3,          // ~17.2 degrees inclination
                                   0.2,          // ~11.5 degrees RAAN
                                   0.1,          // ~5.7 degrees argument of perigee
                                   0.0);         // 0 degrees true anomaly
    
    spdlog::info("\nHyperbolic Orbit Elements:");
    spdlog::info("{}", hyperbolicOrbit.toString());
    
    // Convert to Cartesian
    CartesianState hyperbolicCartState = orbitalElementsToCartesian(hyperbolicOrbit);
    spdlog::info("\nHyperbolic Orbit in Cartesian Coordinates:");
    spdlog::info("{}", hyperbolicCartState.toString());
}

void testOrbitPropagator() {
    spdlog::info("\nTesting Orbit Propagator");
    
    // Create a circular orbit at 7000 km altitude
    Vector3 position(7000000.0 + EARTH_RADIUS, 0.0, 0.0);  // 7000 km altitude along x-axis
    
    // Calculate circular orbit velocity
    double radius = position.norm();
    double circularVelocity = std::sqrt(EARTH_MU / radius);
    Vector3 velocity(0.0, circularVelocity, 0.0);  // Velocity along y-axis for circular orbit
    
    CartesianState initialState(position, velocity);
    spdlog::info("Initial Cartesian State:");
    spdlog::info("{}", initialState.toString());
    
    // Convert to orbital elements to verify
    OrbitalElements initialElements = cartesianToOrbitalElements(initialState);
    spdlog::info("\nInitial Orbital Elements:");
    spdlog::info("{}", initialElements.toString());
    
    // Create propagator
    OrbitPropagator propagator;
    propagator.setInitialState(initialState);
    
    // Set integration parameters
    propagator.setStepSize(60.0);  // 60-second time step
    propagator.setAbsoluteTolerance(1e-10);
    propagator.setRelativeTolerance(1e-8);
    
    // Create storage for propagation history
    std::vector<CartesianState> stateHistory;
    std::vector<double> timeHistory;
    
    // Observer function to record state history
    auto observer = [&stateHistory, &timeHistory](const CartesianState& state, double time) {
        stateHistory.push_back(state);
        timeHistory.push_back(time);
    };
    
    // Propagate for one orbital period
    double orbitalPeriod = initialElements.getPeriod();
    spdlog::info("\nOrbital Period: {:.2f} seconds ({:.2f} minutes)", 
                 orbitalPeriod, orbitalPeriod / 60.0);
    
    // Propagate
    CartesianState finalState = propagator.propagateToTime(orbitalPeriod, observer);
    
    spdlog::info("\nFinal Cartesian State after one orbit:");
    spdlog::info("{}", finalState.toString());
    
    // Convert to orbital elements to verify
    OrbitalElements finalElements = cartesianToOrbitalElements(finalState);
    spdlog::info("\nFinal Orbital Elements:");
    spdlog::info("{}", finalElements.toString());
    
    // Check the difference between initial and final states
    Vector3 positionDifference = finalState.getPosition() - initialState.getPosition();
    Vector3 velocityDifference = finalState.getVelocity() - initialState.getVelocity();
    
    spdlog::info("\nDifference after one orbit:");
    spdlog::info("Position difference: {} m (magnitude: {:.2f} m)", 
                 vectorToString(positionDifference), positionDifference.norm());
    spdlog::info("Velocity difference: {} m/s (magnitude: {:.6f} m/s)", 
                 vectorToString(velocityDifference), velocityDifference.norm());
    
    // Log a few points from the state history
    spdlog::info("\nState history sample points:");
    int numPoints = std::min(5, static_cast<int>(stateHistory.size()));
    for (int i = 0; i < numPoints; ++i) {
        int index = static_cast<int>(stateHistory.size() - 1) * i / (numPoints - 1);
        spdlog::info("Time: {:.2f} s, Position: {}", 
                     timeHistory[index], vectorToString(stateHistory[index].getPosition()));
    }
    
    // Test an elliptical orbit
    spdlog::info("\n\nTesting propagation of an elliptical orbit:");
    OrbitalElements ellipticalElements(8000000.0,  // 8000 km semi-major axis
                                      0.2,        // 0.2 eccentricity
                                      0.0,        // 0 inclination (equatorial)
                                      0.0,        // 0 RAAN
                                      0.0,        // 0 argument of perigee
                                      0.0);       // 0 true anomaly (start at periapsis)
    
    spdlog::info("Elliptical Orbital Elements:");
    spdlog::info("{}", ellipticalElements.toString());
    
    // Convert to Cartesian
    CartesianState ellipticalInitialState = orbitalElementsToCartesian(ellipticalElements);
    spdlog::info("\nElliptical Initial Cartesian State:");
    spdlog::info("{}", ellipticalInitialState.toString());
    
    // Create a new propagator for elliptical orbit
    OrbitPropagator ellipticalPropagator;
    ellipticalPropagator.setInitialState(ellipticalInitialState);
    
    // Clear history
    stateHistory.clear();
    timeHistory.clear();
    
    // Propagate for one orbital period
    double ellipticalPeriod = ellipticalElements.getPeriod();
    spdlog::info("\nElliptical Orbital Period: {:.2f} seconds ({:.2f} minutes)", 
                 ellipticalPeriod, ellipticalPeriod / 60.0);
    
    // Propagate
    CartesianState ellipticalFinalState = ellipticalPropagator.propagateToTime(ellipticalPeriod, observer);
    
    spdlog::info("\nElliptical Final Cartesian State after one orbit:");
    spdlog::info("{}", ellipticalFinalState.toString());
    
    // Convert back to elements
    OrbitalElements ellipticalFinalElements = cartesianToOrbitalElements(ellipticalFinalState);
    spdlog::info("\nElliptical Final Orbital Elements:");
    spdlog::info("{}", ellipticalFinalElements.toString());
    
    // Check the difference between initial and final states
    Vector3 ellipticalPositionDiff = ellipticalFinalState.getPosition() - ellipticalInitialState.getPosition();
    Vector3 ellipticalVelocityDiff = ellipticalFinalState.getVelocity() - ellipticalInitialState.getVelocity();
    
    spdlog::info("\nElliptical Orbit Difference after one orbit:");
    spdlog::info("Position difference: {} m (magnitude: {:.2f} m)", 
                 vectorToString(ellipticalPositionDiff), ellipticalPositionDiff.norm());
    spdlog::info("Velocity difference: {} m/s (magnitude: {:.6f} m/s)", 
                 vectorToString(ellipticalVelocityDiff), ellipticalVelocityDiff.norm());
}

void testPerturbations() {
    spdlog::info("\nTesting Orbital Perturbations");
    
    // Create a circular orbit at 500 km altitude (LEO)
    double altitude = 500000.0;  // 500 km in meters
    Vector3 position(EARTH_RADIUS + altitude, 0.0, 0.0);  // Initial position along x-axis
    
    // Calculate circular orbit velocity
    double radius = position.norm();
    double circularVelocity = std::sqrt(EARTH_MU / radius);
    Vector3 velocity(0.0, circularVelocity, 0.0);  // Velocity along y-axis for circular orbit
    
    CartesianState initialState(position, velocity);
    spdlog::info("Initial Cartesian State (500 km altitude):");
    spdlog::info("{}", initialState.toString());
    
    // Convert to orbital elements to verify
    OrbitalElements initialElements = cartesianToOrbitalElements(initialState);
    spdlog::info("\nInitial Orbital Elements:");
    spdlog::info("{}", initialElements.toString());
    
    // Create a base propagator without perturbations
    OrbitPropagator basePropagator;
    basePropagator.setInitialState(initialState);
    basePropagator.setStepSize(30.0);  // 30-second time step
    
    // Create propagator with J2 perturbation
    OrbitPropagator j2Propagator;
    j2Propagator.setInitialState(initialState);
    j2Propagator.setStepSize(30.0);
    
    // Add J2 perturbation
    auto j2Perturbation = std::make_shared<J2Perturbation>();
    j2Propagator.addPerturbation(j2Perturbation);
    
    // Create propagator with all perturbations
    OrbitPropagator fullPropagator;
    fullPropagator.setInitialState(initialState);
    fullPropagator.setStepSize(30.0);
    
    // Add all perturbations
    fullPropagator.addPerturbation(std::make_shared<J2Perturbation>());
    
    // Create drag perturbation for a typical satellite
    auto drag = std::make_shared<AtmosphericDrag>(
        2.2,       // Cd: 2.2 drag coefficient
        10.0,      // A: 10 m² cross-sectional area
        1000.0     // m: 1000 kg mass
    );
    fullPropagator.addPerturbation(drag);
    
    // Create solar radiation pressure perturbation
    auto srp = std::make_shared<SolarRadiationPressure>(
        1.4,       // Cr: 1.4 reflectivity coefficient
        10.0,      // A: 10 m² cross-sectional area
        1000.0     // m: 1000 kg mass
    );
    fullPropagator.addPerturbation(srp);
    
    // Add third-body perturbations
    fullPropagator.addPerturbation(std::make_shared<ThirdBodyPerturbation>(ThirdBodyPerturbation::Body::SUN));
    fullPropagator.addPerturbation(std::make_shared<ThirdBodyPerturbation>(ThirdBodyPerturbation::Body::MOON));
    
    // Propagation time - 1 day
    double propagationTime = 24.0 * 3600.0; // 24 hours in seconds
    
    // Storage for state history
    struct StateHistory {
        std::vector<double> times;
        std::vector<Vector3> positions;
        std::vector<Vector3> velocities;
    };
    
    StateHistory baseHistory, j2History, fullHistory;
    
    // Observer functions
    auto baseObserver = [&baseHistory](const CartesianState& state, double time) {
        baseHistory.times.push_back(time);
        baseHistory.positions.push_back(state.getPosition());
        baseHistory.velocities.push_back(state.getVelocity());
    };
    
    auto j2Observer = [&j2History](const CartesianState& state, double time) {
        j2History.times.push_back(time);
        j2History.positions.push_back(state.getPosition());
        j2History.velocities.push_back(state.getVelocity());
    };
    
    auto fullObserver = [&fullHistory](const CartesianState& state, double time) {
        fullHistory.times.push_back(time);
        fullHistory.positions.push_back(state.getPosition());
        fullHistory.velocities.push_back(state.getVelocity());
    };
    
    // Propagate
    spdlog::info("\nPropagating orbit for 24 hours...");
    
    CartesianState baseFinalState = basePropagator.propagateToTime(propagationTime, baseObserver);
    CartesianState j2FinalState = j2Propagator.propagateToTime(propagationTime, j2Observer);
    CartesianState fullFinalState = fullPropagator.propagateToTime(propagationTime, fullObserver);
    
    // Convert to orbital elements
    OrbitalElements baseElements = cartesianToOrbitalElements(baseFinalState);
    OrbitalElements j2Elements = cartesianToOrbitalElements(j2FinalState);
    OrbitalElements fullElements = cartesianToOrbitalElements(fullFinalState);
    
    // Output results
    spdlog::info("\nAfter 24 hours of propagation:");
    
    spdlog::info("\nBase (Two-Body) Final Orbital Elements:");
    spdlog::info("{}", baseElements.toString());
    
    spdlog::info("\nJ2 Perturbation Final Orbital Elements:");
    spdlog::info("{}", j2Elements.toString());
    
    spdlog::info("\nFull Perturbation Model Final Orbital Elements:");
    spdlog::info("{}", fullElements.toString());
    
    // Calculate and output the differences
    double semiMajorDiff = j2Elements.getSemiMajorAxis() - baseElements.getSemiMajorAxis();
    double eccentricityDiff = j2Elements.getEccentricity() - baseElements.getEccentricity();
    double inclinationDiff = (j2Elements.getInclination() - baseElements.getInclination()) * 180.0 / M_PI;
    double raanDiff = (j2Elements.getRightAscension() - baseElements.getRightAscension()) * 180.0 / M_PI;
    
    spdlog::info("\nJ2 Perturbation Effects (differences from two-body model):");
    spdlog::info("Semi-Major Axis Change: {:.2f} m", semiMajorDiff);
    spdlog::info("Eccentricity Change: {:.8f}", eccentricityDiff);
    spdlog::info("Inclination Change: {:.6f} deg", inclinationDiff);
    spdlog::info("RAAN Change: {:.6f} deg", raanDiff);
    
    // Save data to CSV for external plotting
    std::ofstream dataFile("orbit_propagation_data.csv");
    dataFile << "Time,X_Base,Y_Base,Z_Base,X_J2,Y_J2,Z_J2,X_Full,Y_Full,Z_Full\n";
    
    // Determine how many data points to output (to keep the file size reasonable)
    int numDataPoints = std::min(static_cast<int>(baseHistory.times.size()), 1000);
    int step = std::max(1, static_cast<int>(baseHistory.times.size()) / numDataPoints);
    
    for (int i = 0; i < baseHistory.times.size(); i += step) {
        dataFile << baseHistory.times[i] << ","
                 << baseHistory.positions[i].x() << ","
                 << baseHistory.positions[i].y() << ","
                 << baseHistory.positions[i].z() << ",";
        
        // Find the closest time point in the J2 history
        int j2Index = i;
        if (j2Index < j2History.times.size()) {
            dataFile << j2History.positions[j2Index].x() << ","
                     << j2History.positions[j2Index].y() << ","
                     << j2History.positions[j2Index].z() << ",";
        } else {
            dataFile << ",,," ;
        }
        
        // Find the closest time point in the full history
        int fullIndex = i;
        if (fullIndex < fullHistory.times.size()) {
            dataFile << fullHistory.positions[fullIndex].x() << ","
                     << fullHistory.positions[fullIndex].y() << ","
                     << fullHistory.positions[fullIndex].z();
        } else {
            dataFile << ",,,";
        }
        
        dataFile << "\n";
    }
    
    dataFile.close();
    spdlog::info("\nSaved propagation data to orbit_propagation_data.csv");
    
    // Test effect of different perturbations separately
    spdlog::info("\nTesting Each Perturbation Separately on Orbital Elements:");
    
    // Create separate propagators for each perturbation
    std::vector<std::pair<std::string, OrbitPropagator>> propagators;
    
    // Two-body only (reference)
    OrbitPropagator twoBodyProp;
    twoBodyProp.setInitialState(initialState);
    twoBodyProp.setStepSize(30.0);
    propagators.push_back({"Two-Body", twoBodyProp});
    
    // J2 only
    OrbitPropagator j2OnlyProp;
    j2OnlyProp.setInitialState(initialState);
    j2OnlyProp.setStepSize(30.0);
    j2OnlyProp.addPerturbation(std::make_shared<J2Perturbation>());
    propagators.push_back({"J2", j2OnlyProp});
    
    // Drag only
    OrbitPropagator dragOnlyProp;
    dragOnlyProp.setInitialState(initialState);
    dragOnlyProp.setStepSize(30.0);
    dragOnlyProp.addPerturbation(std::make_shared<AtmosphericDrag>(2.2, 10.0, 1000.0));
    propagators.push_back({"Drag", dragOnlyProp});
    
    // SRP only
    OrbitPropagator srpOnlyProp;
    srpOnlyProp.setInitialState(initialState);
    srpOnlyProp.setStepSize(30.0);
    srpOnlyProp.addPerturbation(std::make_shared<SolarRadiationPressure>(1.4, 10.0, 1000.0));
    propagators.push_back({"SRP", srpOnlyProp});
    
    // Moon only
    OrbitPropagator moonOnlyProp;
    moonOnlyProp.setInitialState(initialState);
    moonOnlyProp.setStepSize(30.0);
    moonOnlyProp.addPerturbation(std::make_shared<ThirdBodyPerturbation>(ThirdBodyPerturbation::Body::MOON));
    propagators.push_back({"Moon", moonOnlyProp});
    
    // Sun only
    OrbitPropagator sunOnlyProp;
    sunOnlyProp.setInitialState(initialState);
    sunOnlyProp.setStepSize(30.0);
    sunOnlyProp.addPerturbation(std::make_shared<ThirdBodyPerturbation>(ThirdBodyPerturbation::Body::SUN));
    propagators.push_back({"Sun", sunOnlyProp});
    
    // All perturbations
    propagators.push_back({"All", fullPropagator});
    
    // Propagate each for 24 hours and display the results
    for (auto& [name, propagator] : propagators) {
        // Propagate
        CartesianState finalState = propagator.propagateToTime(propagationTime);
        
        // Convert to orbital elements
        OrbitalElements elements = cartesianToOrbitalElements(finalState);
        
        // Output key orbital elements
        spdlog::info("\n{} Model after 24 hours:", name);
        spdlog::info("Semi-Major Axis: {:.1f} km", elements.getSemiMajorAxis() / 1000.0);
        spdlog::info("Eccentricity: {:.8f}", elements.getEccentricity());
        spdlog::info("Inclination: {:.6f} deg", elements.getInclination() * 180.0 / M_PI);
        spdlog::info("RAAN: {:.6f} deg", elements.getRightAscension() * 180.0 / M_PI);
    }
}

void testOrbitManeuvers() {
    spdlog::info("\nTesting Orbit Maneuvers");
    
    // Create a circular orbit at 500 km altitude (LEO)
    double initialAltitude = 500000.0;  // 500 km in meters
    double initialRadius = EARTH_RADIUS + initialAltitude;
    Vector3 position(initialRadius, 0.0, 0.0);  // Initial position along x-axis
    
    // Calculate circular orbit velocity
    double circularVelocity = std::sqrt(EARTH_MU / initialRadius);
    Vector3 velocity(0.0, circularVelocity, 0.0);  // Velocity along y-axis for circular orbit
    
    CartesianState initialState(position, velocity);
    spdlog::info("Initial Cartesian State (500 km altitude):");
    spdlog::info("{}", initialState.toString());
    
    // Convert to orbital elements to verify
    OrbitalElements initialElements = cartesianToOrbitalElements(initialState);
    spdlog::info("\nInitial Orbital Elements:");
    spdlog::info("{}", initialElements.toString());
    
    // Test 1: Hohmann Transfer to Higher Orbit
    spdlog::info("\nTest 1: Hohmann Transfer from 500 km to 1000 km");
    
    double finalAltitude = 1000000.0;  // 1000 km in meters
    double finalRadius = EARTH_RADIUS + finalAltitude;
    
    // Calculate delta-V for Hohmann transfer
    double hohmannDeltaV = calculateHohmannTransferDeltaV(initialRadius, finalRadius);
    double hohmannTimeOfFlight = calculateHohmannTimeOfFlight(initialRadius, finalRadius);
    
    spdlog::info("Hohmann Transfer Delta-V: {:.2f} m/s", hohmannDeltaV);
    spdlog::info("Hohmann Transfer Time of Flight: {:.2f} seconds ({:.2f} minutes)",
                hohmannTimeOfFlight, hohmannTimeOfFlight / 60.0);
    
    // Create maneuvers for Hohmann transfer
    auto initialBurn = createHohmannTransferInitialBurn(initialRadius, finalRadius);
    auto finalBurn = createHohmannTransferFinalBurn(initialRadius, finalRadius, hohmannTimeOfFlight);
    
    spdlog::info("Initial Burn Delta-V: {:.2f} m/s", initialBurn->getDeltaV());
    spdlog::info("Final Burn Delta-V: {:.2f} m/s", finalBurn->getDeltaV());
    
    // Create propagator
    OrbitPropagator hohmannPropagator;
    hohmannPropagator.setInitialState(initialState);
    hohmannPropagator.setStepSize(30.0);  // 30-second time step
    
    // Add maneuvers
    hohmannPropagator.addManeuver(initialBurn);
    hohmannPropagator.addManeuver(finalBurn);
    
    // Storage for propagation history
    std::vector<CartesianState> hohmannStateHistory;
    std::vector<double> hohmannTimeHistory;
    
    // Observer function
    auto observer = [&hohmannStateHistory, &hohmannTimeHistory](const CartesianState& state, double time) {
        hohmannStateHistory.push_back(state);
        hohmannTimeHistory.push_back(time);
    };
    
    // Propagate for the Hohmann transfer plus half an orbit
    double propagationTime = hohmannTimeOfFlight + 0.5 * 2.0 * M_PI * std::sqrt(std::pow(finalRadius, 3) / EARTH_MU);
    CartesianState finalState = hohmannPropagator.propagateToTime(propagationTime, observer);
    
    // Convert to orbital elements
    OrbitalElements finalElements = cartesianToOrbitalElements(finalState);
    
    spdlog::info("\nFinal Orbital Elements after Hohmann Transfer:");
    spdlog::info("{}", finalElements.toString());
    spdlog::info("Final Altitude: {:.2f} km", (finalElements.getSemiMajorAxis() - EARTH_RADIUS) / 1000.0);
    
    // Save Hohmann transfer data for plotting
    std::ofstream hohmannDataFile("hohmann_transfer_data.csv");
    hohmannDataFile << "Time,X,Y,Z,Radius\n";
    
    for (size_t i = 0; i < hohmannStateHistory.size(); ++i) {
        const auto& state = hohmannStateHistory[i];
        double time = hohmannTimeHistory[i];
        double radius = state.getPosition().norm();
        
        hohmannDataFile << time << ","
                       << state.getPosition().x() << ","
                       << state.getPosition().y() << ","
                       << state.getPosition().z() << ","
                       << radius << "\n";
    }
    
    hohmannDataFile.close();
    spdlog::info("Saved Hohmann transfer data to hohmann_transfer_data.csv");
    
    // Test 2: Plane Change Maneuver
    spdlog::info("\n\nTest 2: Plane Change Maneuver");
    
    // Create a circular orbit in the equatorial plane
    CartesianState equatorialState(position, velocity);
    
    // Create propagator
    OrbitPropagator planeChangePropagator;
    planeChangePropagator.setInitialState(equatorialState);
    planeChangePropagator.setStepSize(30.0);
    
    // Define plane change angle (30 degrees)
    double planeChangeAngle = 30.0 * M_PI / 180.0;  // radians
    
    // Calculate delta-V for plane change
    double planeChangeDeltaV = calculatePlaneChangeDeltaV(circularVelocity, planeChangeAngle);
    spdlog::info("Plane Change Delta-V for {:.1f} degrees: {:.2f} m/s", 
                planeChangeAngle * 180.0 / M_PI, planeChangeDeltaV);
    
    // Create plane change maneuver
    // Execute after half an orbit to ensure we're at the ascending/descending node
    double halfOrbitTime = M_PI * std::sqrt(std::pow(initialRadius, 3) / EARTH_MU);
    auto planeChangeManeuver = createPlaneChangeManeuver(circularVelocity, planeChangeAngle, halfOrbitTime);
    
    // Add maneuver
    planeChangePropagator.addManeuver(planeChangeManeuver);
    
    // Storage for propagation history
    std::vector<CartesianState> planeChangeStateHistory;
    std::vector<double> planeChangeTimeHistory;
    
    // Observer function
    auto planeChangeObserver = [&planeChangeStateHistory, &planeChangeTimeHistory](const CartesianState& state, double time) {
        planeChangeStateHistory.push_back(state);
        planeChangeTimeHistory.push_back(time);
    };
    
    // Propagate for one full orbit after the plane change
    double fullOrbitTime = 2.0 * M_PI * std::sqrt(std::pow(initialRadius, 3) / EARTH_MU);
    CartesianState planeChangeFinalState = planeChangePropagator.propagateToTime(halfOrbitTime + fullOrbitTime, planeChangeObserver);
    
    // Convert to orbital elements
    OrbitalElements planeChangeFinalElements = cartesianToOrbitalElements(planeChangeFinalState);
    
    spdlog::info("\nFinal Orbital Elements after Plane Change:");
    spdlog::info("{}", planeChangeFinalElements.toString());
    spdlog::info("Final Inclination: {:.2f} degrees", planeChangeFinalElements.getInclination() * 180.0 / M_PI);
    
    // Save plane change data for plotting
    std::ofstream planeChangeDataFile("plane_change_data.csv");
    planeChangeDataFile << "Time,X,Y,Z,Inclination\n";
    
    for (size_t i = 0; i < planeChangeStateHistory.size(); ++i) {
        const auto& state = planeChangeStateHistory[i];
        double time = planeChangeTimeHistory[i];
        
        // Calculate inclination at each step
        OrbitalElements elements = cartesianToOrbitalElements(state);
        double inclination = elements.getInclination() * 180.0 / M_PI;
        
        planeChangeDataFile << time << ","
                           << state.getPosition().x() << ","
                           << state.getPosition().y() << ","
                           << state.getPosition().z() << ","
                           << inclination << "\n";
    }
    
    planeChangeDataFile.close();
    spdlog::info("Saved plane change data to plane_change_data.csv");
    
    // Test 3: Finite Burn Maneuver
    spdlog::info("\n\nTest 3: Finite Burn Maneuver");
    
    // Create a new state for testing finite burns
    CartesianState finiteBurnState(position, velocity);
    
    // Create propagator
    OrbitPropagator finiteBurnPropagator;
    finiteBurnPropagator.setInitialState(finiteBurnState);
    finiteBurnPropagator.setStepSize(1.0);  // 1-second time step for better resolution
    
    // Define satellite parameters
    double satelliteMass = 1000.0;      // 1000 kg
    double thrust = 10.0;               // 10 N
    double specificImpulse = 300.0;     // 300 s (typical for chemical propulsion)
    
    // Create a prograde finite burn maneuver for 500 seconds
    double burnDuration = 500.0;  // seconds
    auto finiteBurnManeuver = std::make_shared<FiniteBurnManeuver>(
        thrust, specificImpulse, satelliteMass, burnDuration, 
        OrbitManeuver::Direction::PROGRADE, 0.0);
    
    // Calculate the expected delta-V
    double expectedDeltaV = finiteBurnManeuver->calculateTotalDeltaV();
    spdlog::info("Finite Burn Expected Delta-V: {:.2f} m/s", expectedDeltaV);
    spdlog::info("Finite Burn Duration: {:.1f} seconds", burnDuration);
    spdlog::info("Thrust: {:.1f} N, Specific Impulse: {:.1f} s", thrust, specificImpulse);
    spdlog::info("Initial Mass: {:.1f} kg", satelliteMass);
    
    // Calculate the final mass after the burn
    double finalMass = finiteBurnManeuver->calculateCurrentMass(burnDuration);
    spdlog::info("Fuel Consumed: {:.2f} kg", satelliteMass - finalMass);
    spdlog::info("Final Mass: {:.2f} kg", finalMass);
    
    // Add maneuver
    finiteBurnPropagator.addManeuver(finiteBurnManeuver);
    
    // Storage for propagation history
    std::vector<CartesianState> finiteBurnStateHistory;
    std::vector<double> finiteBurnTimeHistory;
    
    // Observer function
    auto finiteBurnObserver = [&finiteBurnStateHistory, &finiteBurnTimeHistory](const CartesianState& state, double time) {
        finiteBurnStateHistory.push_back(state);
        finiteBurnTimeHistory.push_back(time);
    };
    
    // Propagate for the burn duration plus half an orbit
    double halfOrbitAfterBurn = burnDuration + 0.5 * fullOrbitTime;
    CartesianState finiteBurnFinalState = finiteBurnPropagator.propagateToTime(halfOrbitAfterBurn, finiteBurnObserver);
    
    // Convert to orbital elements
    OrbitalElements finiteBurnFinalElements = cartesianToOrbitalElements(finiteBurnFinalState);
    
    spdlog::info("\nFinal Orbital Elements after Finite Burn:");
    spdlog::info("{}", finiteBurnFinalElements.toString());
    spdlog::info("Final Semi-Major Axis: {:.2f} km", finiteBurnFinalElements.getSemiMajorAxis() / 1000.0);
    spdlog::info("Final Altitude: {:.2f} km", (finiteBurnFinalElements.getSemiMajorAxis() - EARTH_RADIUS) / 1000.0);
    
    // Calculate the actual delta-V achieved
    double initialVelocity = velocity.norm();
    double finalVelocity = finiteBurnFinalState.getVelocity().norm();
    
    spdlog::info("Initial Velocity: {:.2f} m/s", initialVelocity);
    spdlog::info("Final Velocity: {:.2f} m/s", finalVelocity);
    spdlog::info("Velocity Change: {:.2f} m/s", finalVelocity - initialVelocity);
    
    // Save finite burn data for plotting
    std::ofstream finiteBurnDataFile("finite_burn_data.csv");
    finiteBurnDataFile << "Time,X,Y,Z,Velocity,Altitude\n";
    
    for (size_t i = 0; i < finiteBurnStateHistory.size(); ++i) {
        const auto& state = finiteBurnStateHistory[i];
        double time = finiteBurnTimeHistory[i];
        double velocity = state.getVelocity().norm();
        double altitude = (state.getPosition().norm() - EARTH_RADIUS);
        
        finiteBurnDataFile << time << ","
                          << state.getPosition().x() << ","
                          << state.getPosition().y() << ","
                          << state.getPosition().z() << ","
                          << velocity << ","
                          << altitude << "\n";
    }
    
    finiteBurnDataFile.close();
    spdlog::info("Saved finite burn data to finite_burn_data.csv");
    
    // Test 4: Combined Maneuver (Hohmann Transfer + Plane Change)
    spdlog::info("\n\nTest 4: Combined Maneuver (Hohmann Transfer + Plane Change)");
    
    // Calculate delta-V for combined maneuver
    double combinedDeltaV = calculateCombinedManeuverDeltaV(initialRadius, finalRadius, planeChangeAngle);
    spdlog::info("Combined Maneuver Delta-V: {:.2f} m/s", combinedDeltaV);
    spdlog::info("Separate Maneuvers Delta-V: {:.2f} m/s", hohmannDeltaV + planeChangeDeltaV);
    spdlog::info("Delta-V Savings: {:.2f} m/s", hohmannDeltaV + planeChangeDeltaV - combinedDeltaV);
}

int main() {
    spdlog::info("Starting Satellite Simulator Physics Engine");
    
    // Test vector math utilities
    testVectorMath();
    
    // Test orbital state representations
    testOrbitalState();
    
    // Test basic orbit propagator
    testOrbitPropagator();
    
    // Test perturbations
    testPerturbations();
    
    // Test orbit maneuvers
    testOrbitManeuvers();
    
    spdlog::info("Physics Engine test completed successfully");
    return 0;
}