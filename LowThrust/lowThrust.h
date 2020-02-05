/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/getRecommendedBaseFunctionsHodographicShaping.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "../applicationOutput.h"

using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::aerodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::mathematical_constants;
using namespace tudat::reference_frames;
using namespace tudat::shape_based_methods;
using namespace tudat::low_thrust_trajectories;
using namespace tudat;

namespace tudat_applications
{
namespace PropagationOptimization2020
{

double getTrajectoryTimeOfFlight( const std::vector< double >& trajectoryParameters);

double getTrajectoryInitialTime( const std::vector< double >& trajectoryParameters, const double bufferTime = 0.0 );

double getTrajectoryFinalTime( const std::vector< double >& trajectoryParameters, const double bufferTime = 0.0 );


//! Get the propagation termination settings for the lunar ascent.
/*!
 * \param simulationStartEpoch Start time of the simulation in seconds.
 * \param maximumDuration Time in seconds, specifying the maximum time duration before which the
 * simulation should stop.
 * \param terminationAltitude Altitude in meters, specifying the maximum altitude before which the
 * simulation should stop.
 * \param vehicleDryMass Dry mass of the vehicle in kg. This is value is used to create a termination
 * condition that mandates the simulation to stop once all fuel has been used up.
 * \return Shared pointer to the PropagationTerminationSettings object.
 */
std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings(
        const std::vector< double >& trajectoryParameters,
        const double targetDistance,
        const double trajectoryTimeBuffer );

std::shared_ptr< HodographicShaping > createHodographicShapingObject(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap );

std::shared_ptr< ThrustAccelerationSettings > getThrustAccelerationSettingsFromParameters(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap );

Eigen::Vector6d getHodographicLowThrustStateAtEpoch(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap,
        const double evaluationTime );

class LowThrustProblem
{
public:

    // Constructor that we will actually use
    LowThrustProblem(
            const simulation_setup::NamedBodyMap bodyMap,
            const std::shared_ptr< IntegratorSettings< > > integratorSettings,
            const std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings,
            double specificImpulse,
            double minimumMarsDistance,
            double timeBuffer,
            const bool performPropagation = true );

    // Standard constructor
    LowThrustProblem( )
    {
    }

    //! Function to retrieve the map with the propagated state history of the last run
    std::map< double, Eigen::VectorXd > getLastRunPropagatedStateHistory( ) const
    {
        return dynamicsSimulator_->getEquationsOfMotionNumericalSolution( );
    }

    //! Function to retrieve the map with the dependent variable history of the last run
    std::map< double, Eigen::VectorXd > getLastRunDependentVariableHistory( ) const
    {
        return dynamicsSimulator_->getDependentVariableHistory( );
    }

    //! Function to retrieve a shared pointer to the dynamics simulator of the last run
    std::shared_ptr< SingleArcDynamicsSimulator< > > getLastRunDynamicsSimulator( )
    {
        return dynamicsSimulator_;
    }


    //! Fitness function, called to run the simulation. In this form, it is compatible with
    //! the Pagmo optimization library.
    //! \param shapeParameters Decision vector, containing the shape parameters for which
    //! the propagation needs to be run, to find the corresponding fitness value.
    //! \return Returns the vector with doubles, describing the fitness belonging
    //! to the shape parameters passed to the function. Currently returns an empty
    //! vector.
    //!
    std::vector< double > fitness( std::vector< double >& x ) const;

    std::shared_ptr< HodographicShaping > getHodographicShaping( )
    {
        return hodographicShaping_;
    }

private:

    //! Instance variable holding the body map for the simulation
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Object holding the integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings_;

    //! Object holding the propagator settings
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings_;

    //! Object holding the translational state propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings_;

    double specificImpulse_;

    double minimumMarsDistance_;

    double timeBuffer_;

    bool performPropagation_;

    mutable std::shared_ptr< HodographicShaping > hodographicShaping_;

    //! Object holding the dynamics simulator
    mutable std::shared_ptr<SingleArcDynamicsSimulator< > > dynamicsSimulator_;

};

} // Namespace tudat_applications
} // Namespace PropagationOptimization2020
