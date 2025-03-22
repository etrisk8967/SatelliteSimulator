// physics_engine/include/perturbation_force.h
#pragma once

#include "vector_math.h"
#include "orbital_state.h"

namespace SatelliteSimulator {

    // Base class for all perturbation forces
    class PerturbationForce {
    public:
        virtual ~PerturbationForce() = default;
        
        // Calculate the acceleration due to this perturbation
        virtual Vector3 calculateAcceleration(const CartesianState& state, double time) const = 0;
        
        // Get name of the perturbation (for logging and debugging)
        virtual std::string getName() const = 0;
        
        // Enable/disable this perturbation
        void setEnabled(bool enabled) { enabled_ = enabled; }
        bool isEnabled() const { return enabled_; }
        
    protected:
        bool enabled_ = true;  // Perturbations are enabled by default
    };

    // J2 Perturbation (Earth's oblateness)
    class J2Perturbation : public PerturbationForce {
    public:
        J2Perturbation(double j2Value = 0.00108263);  // Default J2 value for Earth
        
        Vector3 calculateAcceleration(const CartesianState& state, double time) const override;
        std::string getName() const override { return "J2 Perturbation"; }
        
        // Get/set J2 coefficient
        double getJ2Value() const { return j2Value_; }
        void setJ2Value(double j2Value) { j2Value_ = j2Value; }
        
    private:
        double j2Value_;  // J2 coefficient
    };

    // Atmospheric Drag Perturbation
    class AtmosphericDrag : public PerturbationForce {
    public:
        AtmosphericDrag(double dragCoefficient = 2.2, double crossSectionalArea = 1.0, 
                        double satelliteMass = 100.0);
        
        Vector3 calculateAcceleration(const CartesianState& state, double time) const override;
        std::string getName() const override { return "Atmospheric Drag"; }
        
        // Get/set parameters
        double getDragCoefficient() const { return dragCoefficient_; }
        void setDragCoefficient(double cd) { dragCoefficient_ = cd; }
        
        double getCrossSectionalArea() const { return crossSectionalArea_; }
        void setCrossSectionalArea(double area) { crossSectionalArea_ = area; }
        
        double getSatelliteMass() const { return satelliteMass_; }
        void setSatelliteMass(double mass) { satelliteMass_ = mass; }
        
    private:
        double dragCoefficient_;      // Cd (dimensionless)
        double crossSectionalArea_;   // A (m²)
        double satelliteMass_;        // m (kg)
        
        // Calculate atmospheric density using exponential model
        double calculateAtmosphericDensity(double altitude) const;
    };

    // Solar Radiation Pressure Perturbation
    class SolarRadiationPressure : public PerturbationForce {
    public:
        SolarRadiationPressure(double reflectivityCoefficient = 1.0, 
                              double crossSectionalArea = 1.0,
                              double satelliteMass = 100.0);
        
        Vector3 calculateAcceleration(const CartesianState& state, double time) const override;
        std::string getName() const override { return "Solar Radiation Pressure"; }
        
        // Get/set parameters
        double getReflectivityCoefficient() const { return reflectivityCoefficient_; }
        void setReflectivityCoefficient(double cr) { reflectivityCoefficient_ = cr; }
        
        double getCrossSectionalArea() const { return crossSectionalArea_; }
        void setCrossSectionalArea(double area) { crossSectionalArea_ = area; }
        
        double getSatelliteMass() const { return satelliteMass_; }
        void setSatelliteMass(double mass) { satelliteMass_ = mass; }
        
        // Moved from private to public so it can be accessed by ThirdBodyPerturbation
        Vector3 calculateSunPosition(double time) const;
        
    private:
        double reflectivityCoefficient_;  // Cr (dimensionless)
        double crossSectionalArea_;       // A (m²)
        double satelliteMass_;            // m (kg)
        
        // Check if satellite is in Earth's shadow
        bool isInShadow(const Vector3& satellitePosition, const Vector3& sunPosition) const;
    };

    // Third-Body Perturbation (Moon and Sun)
    class ThirdBodyPerturbation : public PerturbationForce {
    public:
        enum class Body {
            MOON,
            SUN
        };
        
        ThirdBodyPerturbation(Body body);
        
        Vector3 calculateAcceleration(const CartesianState& state, double time) const override;
        std::string getName() const override;
        
        // Get body type
        Body getBodyType() const { return bodyType_; }
        
    private:
        Body bodyType_;
        double bodyMass_;  // Mass of the third body (kg)
        
        // Calculate position of the third body at a given time (simplified model)
        Vector3 calculateBodyPosition(double time) const;
    };

} // namespace SatelliteSimulator