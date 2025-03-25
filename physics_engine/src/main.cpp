#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <memory>
#include "../include/vector_math.h"
#include "../include/orbital_state.h"
#include "../include/orbit_propagator.h"
#include "../include/perturbation_force.h"
using namespace SatelliteSimulator;

struct GeodeticCoordinates {
    double latitude, longitude, altitude; // radians, radians, meters
};

GeodeticCoordinates cartesianToGeodetic(const Vector3& position) {
    double r = position.norm();
    return {
        std::asin(position.z() / r),
        std::atan2(position.y(), position.x()),
        r - EARTH_RADIUS
    };
}

int main() {
    std::cout << "Welcome, space overlord! Let’s launch a satellite!\n";

    // Input Type
    int inputType;
    std::cout << "Input method: (1) Cartesian or (2) Orbital Elements? ";
    std::cin >> inputType;

    CartesianState initialState;
    if (inputType == 1) {
        double x, y, z, vx, vy, vz;
        std::cout << "Position (x, y, z) in meters: ";
        std::cin >> x >> y >> z;
        std::cout << "Velocity (vx, vy, vz) in m/s: ";
        std::cin >> vx >> vy >> vz;
        initialState = CartesianState(Vector3(x, y, z), Vector3(vx, vy, vz));
    } else {
        double a, e, i, raan, omega, nu;
        std::cout << "Semi-major axis (m): "; std::cin >> a;
        std::cout << "Eccentricity: "; std::cin >> e;
        std::cout << "Inclination (deg): "; std::cin >> i;
        std::cout << "RAAN (deg): "; std::cin >> raan;
        std::cout << "Arg of Perigee (deg): "; std::cin >> omega;
        std::cout << "True Anomaly (deg): "; std::cin >> nu;
        OrbitalElements elements(a, e, i * M_PI / 180, raan * M_PI / 180, omega * M_PI / 180, nu * M_PI / 180);
        initialState = orbitalElementsToCartesian(elements);
    }

    // Satellite Parameters
    double Cd, A, m, Cr;
    std::cout << "Drag Coefficient: "; std::cin >> Cd;
    std::cout << "Cross-sectional Area (m²): "; std::cin >> A;
    std::cout << "Mass (kg): "; std::cin >> m;
    std::cout << "Reflectivity Coefficient: "; std::cin >> Cr;

    // Perturbations
    char enable;
    OrbitPropagator prop;
    prop.setInitialState(initialState);
    prop.setStepSize(60.0);
    std::cout << "Enable J2? (Y/N): "; std::cin >> enable;
    if (enable == 'Y') prop.addPerturbation(std::make_shared<J2Perturbation>());
    std::cout << "Enable Drag? (Y/N): "; std::cin >> enable;
    if (enable == 'Y') prop.addPerturbation(std::make_shared<AtmosphericDrag>(Cd, A, m));
    std::cout << "Enable SRP? (Y/N): "; std::cin >> enable;
    if (enable == 'Y') prop.addPerturbation(std::make_shared<SolarRadiationPressure>(Cr, A, m));
    std::cout << "Enable Sun? (Y/N): "; std::cin >> enable;
    if (enable == 'Y') prop.addPerturbation(std::make_shared<ThirdBodyPerturbation>(ThirdBodyPerturbation::Body::SUN));
    std::cout << "Enable Moon? (Y/N): "; std::cin >> enable;
    if (enable == 'Y') prop.addPerturbation(std::make_shared<ThirdBodyPerturbation>(ThirdBodyPerturbation::Body::MOON));

    // Propagation
    bool impacted = false;
    GeodeticCoordinates impactLocation;
    std::vector<CartesianState> states;
    std::vector<double> times;

    // Simplified RK4 propagation
    double t = 0, dt = 60, maxT = 24 * 3600;
    states.push_back(initialState);
    times.push_back(t);
    while (t < maxT && !impacted) {
        CartesianState state = states.back();
        Vector3 k1v = prop.calculateTotalAcceleration(state, t);
        Vector3 k1p = state.getVelocity();

        Vector3 pos2 = state.getPosition() + k1p * dt / 2;
        Vector3 vel2 = state.getVelocity() + k1v * dt / 2;
        Vector3 k2v = prop.calculateTotalAcceleration(CartesianState(pos2, vel2), t + dt / 2);
        Vector3 k2p = vel2;

        Vector3 pos3 = state.getPosition() + k2p * dt / 2;
        Vector3 vel3 = state.getVelocity() + k2v * dt / 2;
        Vector3 k3v = prop.calculateTotalAcceleration(CartesianState(pos3, vel3), t + dt / 2);
        Vector3 k3p = vel3;

        Vector3 pos4 = state.getPosition() + k3p * dt;
        Vector3 vel4 = state.getVelocity() + k3v * dt;
        Vector3 k4v = prop.calculateTotalAcceleration(CartesianState(pos4, vel4), t + dt);
        Vector3 k4p = vel4;

        Vector3 newPos = state.getPosition() + (k1p + 2 * k2p + 2 * k3p + k4p) * dt / 6;
        Vector3 newVel = state.getVelocity() + (k1v + 2 * k2v + 2 * k3v + k4v) * dt / 6;
        CartesianState newState(newPos, newVel);

        if (newPos.norm() <= EARTH_RADIUS) {
            impacted = true;
            impactLocation = cartesianToGeodetic(newPos);
            break;
        }

        states.push_back(newState);
        t += dt;
        times.push_back(t);
    }

    // Output
    std::cout << "Orbit Propagation Results—Hold onto Your Votes!\n";
    std::cout << "Time (s) | Position (m) | Velocity (m/s) | Alt (m) | Lat (deg) | Lon (deg)\n";
    std::cout << "---------|--------------|----------------|---------|-----------|-----------\n";
    for (size_t i = 0; i < states.size(); ++i) {
        auto geo = cartesianToGeodetic(states[i].getPosition());
        std::cout << std::fixed << std::setprecision(0)
                  << times[i] << " | "
                  << vectorToString(states[i].getPosition()) << " | "
                  << vectorToString(states[i].getVelocity()) << " | "
                  << geo.altitude << " | "
                  << std::setprecision(2) << geo.latitude * 180 / M_PI << " | "
                  << geo.longitude * 180 / M_PI << "\n";
    }

    if (impacted) {
        std::cout << "IMPACT ALERT! Time: " << t << " s, Lat: " << impactLocation.latitude * 180 / M_PI
                  << "°, Lon: " << impactLocation.longitude * 180 / M_PI << "°—Send the cleanup crew!\n";
    } else {
        OrbitalElements finalElements = cartesianToOrbitalElements(states.back());
        std::cout << "\nFinal Orbital Elements:\n" << finalElements.toString() << "\n";
        std::cout << "Total Propagation Time: " << maxT << " s\n";
    }

    std::cout << "Mission complete—unlike Congress, we got something done!\n";
    return 0;
}
