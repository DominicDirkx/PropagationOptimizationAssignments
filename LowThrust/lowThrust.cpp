/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "lowThrust.h"

namespace tudat_applications
{
namespace PropagationOptimization2020
{

double getTrajectoryTimeOfFlight( const std::vector< double >& trajectoryParameters)
{
    return trajectoryParameters.at( 1 ) * physical_constants::JULIAN_DAY;
}

double getTrajectoryInitialTime( const std::vector< double >& trajectoryParameters, const double bufferTime  )
{
    return trajectoryParameters.at( 0 ) * physical_constants::JULIAN_DAY + bufferTime;
}

double getTrajectoryFinalTime( const std::vector< double >& trajectoryParameters, const double bufferTime )
{
    return ( trajectoryParameters.at( 0 ) + trajectoryParameters.at( 1 )  )* physical_constants::JULIAN_DAY - bufferTime;
}


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
        const double trajectoryTimeBuffer )
{

    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           getTrajectoryFinalTime( trajectoryParameters, trajectoryTimeBuffer ) ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               relative_distance_dependent_variable, "Vehicle", "Mars" ), targetDistance, true ) );
    return std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

}

std::shared_ptr< HodographicShaping > createHodographicShapingObject(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap )
{
    double initialTime = getTrajectoryInitialTime( trajectoryParameters );
    double timeOfFlight = getTrajectoryTimeOfFlight( trajectoryParameters );
    double finalTime = getTrajectoryFinalTime( trajectoryParameters );

    int numberOfRevolutions = trajectoryParameters.at( 2 );

    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Define settings for the two additional base functions for the radial velocity composite function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    shape_based_methods::getRecommendedRadialVelocityBaseFunctions(
                radialVelocityFunctionComponents, freeCoefficientsRadialVelocityFunction, timeOfFlight );

    {
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );
        freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
        freeCoefficientsRadialVelocityFunction( 0 ) = trajectoryParameters.at( 3 );
        freeCoefficientsRadialVelocityFunction( 1 ) = trajectoryParameters.at( 4 );
    }

    // Get recommended base functions for the normal velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    shape_based_methods::getRecommendedNormalAxialBaseFunctions( normalVelocityFunctionComponents, freeCoefficientsNormalVelocityFunction, timeOfFlight );

    {
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthNormalVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthNormalVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );
        normalVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fourthNormalVelocityBaseFunctionSettings ) );
        normalVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fifthNormalVelocityBaseFunctionSettings ) );
        freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 2 );
        freeCoefficientsNormalVelocityFunction( 0 ) = trajectoryParameters.at( 5 );
        freeCoefficientsNormalVelocityFunction( 1 ) = trajectoryParameters.at( 6 );
    }

    // Get recommended base functions for the axial velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    shape_based_methods::getRecommendedAxialVelocityBaseFunctions(
                axialVelocityFunctionComponents, freeCoefficientsAxialVelocityFunction,
                timeOfFlight, numberOfRevolutions );

    {
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >
                ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fourthAxialVelocityBaseFunctionSettings ) );
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fifthAxialVelocityBaseFunctionSettings ) );

        freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 2 );
        freeCoefficientsAxialVelocityFunction( 0 ) = trajectoryParameters.at( 7 );
        freeCoefficientsAxialVelocityFunction( 1 ) =  trajectoryParameters.at( 8 );
    }

    Eigen::Vector6d initialState = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( initialTime );
    Eigen::Vector6d finalState = bodyMap.at( "Mars" )->getStateInBaseFrameFromEphemeris( finalTime );
    double centralBodyGravitationalParameter = bodyMap.at( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( );

    return std::make_shared< HodographicShaping >(
                initialState, finalState, timeOfFlight, centralBodyGravitationalParameter, numberOfRevolutions,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );
}

std::shared_ptr< ThrustAccelerationSettings > getThrustAccelerationSettingsFromParameters(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap )
{
    return createHodographicShapingObject( trajectoryParameters, bodyMap )->getLowThrustAccelerationSettings(
                bodyMap, "Vehicle", nullptr, nullptr, getTrajectoryInitialTime( trajectoryParameters ) );
}

Eigen::Vector6d getHodographicLowThrustStateAtEpoch(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap,
        const double evaluationTime )
{
    std::map< double, Eigen::Vector6d > stateMap;
    std::vector< double > epochs = { evaluationTime - getTrajectoryInitialTime( trajectoryParameters, 0.0 ) };
    createHodographicShapingObject( trajectoryParameters, bodyMap )->getTrajectory( epochs, stateMap );
    return stateMap.begin( )->second;

}

}

}

using namespace tudat_applications::PropagationOptimization2020;

LowThrustProblem::LowThrustProblem( const simulation_setup::NamedBodyMap bodyMap,
                                    const std::shared_ptr< IntegratorSettings< > > integratorSettings,
                                    const std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings,
                                    double specificImpulse,
                                    double minimumMarsDistance,
                                    double timeBuffer ):
    bodyMap_(bodyMap), integratorSettings_(integratorSettings), propagatorSettings_(propagatorSettings),
    specificImpulse_( specificImpulse ), minimumMarsDistance_( minimumMarsDistance ), timeBuffer_( timeBuffer )
{
    translationalStatePropagatorSettings_ =
            std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< double > >(
                propagatorSettings_->propagatorSettingsMap_.at( translational_state ).at( 0 ) );
}

std::vector< double > LowThrustProblem::fitness( std::vector< double >& trajectoryParameters ) const
{

    double initialPropagationTime = getTrajectoryInitialTime( trajectoryParameters, timeBuffer_ );
    integratorSettings_->initialTime_ = initialPropagationTime;

    // Create trajectory shape object
    hodographicShaping_ = createHodographicShapingObject(
                trajectoryParameters, bodyMap_ );

    // Extract existing acceleration settings, and clear existing self-exerted accelerations of vehicle
    simulation_setup::SelectedAccelerationMap accelerationSettings =
            translationalStatePropagatorSettings_->getAccelerationSettingsMap( );
    accelerationSettings[ "Vehicle" ][ "Vehicle" ].clear( );

    // Retrieve new acceleration model for thrust and set in list of settings
    std::shared_ptr< AccelerationSettings > newThrustSettings = hodographicShaping_->getLowThrustAccelerationSettings(
                bodyMap_, "Vehicle", [=](const double){return specificImpulse_;}, integratorSettings_, getTrajectoryInitialTime( trajectoryParameters ) );
    accelerationSettings[ "Vehicle" ][ "Vehicle" ].push_back( newThrustSettings );

    // Update translational propagatot settings
    translationalStatePropagatorSettings_->resetAccelerationModelsMap(
                accelerationSettings, bodyMap_ );
    Eigen::Vector6d systemInitialState = getHodographicLowThrustStateAtEpoch(
                trajectoryParameters, bodyMap_, initialPropagationTime );
    translationalStatePropagatorSettings_->resetInitialStates( systemInitialState );

    // Update full propagator settings
    propagatorSettings_->resetIntegratedStateModels( bodyMap_ );
    propagatorSettings_->resetInitialStates(
                createCombinedInitialState< double >( propagatorSettings_->propagatorSettingsMap_ ).segment( 0, 7 ) );
    propagatorSettings_->resetTerminationSettings(
                getPropagationTerminationSettings(
                    trajectoryParameters, minimumMarsDistance_, 0.0 ) );

    dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< > >( bodyMap_, integratorSettings_, propagatorSettings_ );

    return {0.0};

}
