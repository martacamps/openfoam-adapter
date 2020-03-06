#ifndef FF_H
#define FF_H

#include "Interface.H"

#include "FF/Velocity.H"
#include "FF/ReDiag.H"
#include "FF/Pressure.H"
#include "FF/Epsilon.H"
#include "FF/PressureGradient.H"
#include "FF/VelocityGradient.H"

#include "fvCFD.H"

// YAML reader - Used to read the adapter's configuration file.
#include "yaml-cpp/yaml.h"

namespace preciceAdapter
{
namespace FF
{

class FluidFluid
{

protected:

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    // TODO: Create a struct for all the parameter names

    //- Solver type
    std::string solverType_ = "none";

    //- Name of the velocity field
    std::string nameU_ = "U";
    
    //- Name of the pressure field
    std::string nameP_ = "p";
    
    //- Name of the turbulent dissipation field
    std::string nameE_ = "epsilon";
    
    //- Name of the Reynolds stress field (the simulations needs to use the turbulenceFields function object)
    std::string nameR_ = "turbulenceProperties:R";

    //- Determine the solver type
    std::string determineSolverType();

    //- Read the FF-related options from the adapter's configuration file
    bool readConfig(const YAML::Node adapterConfig);
    
    //- Inlet volumetric flow rate. The value by default deactivates the mass correction
    // TODO: Should I change this so that it could work with negative flow rates?
    double vDot_ = -1;

public:

    //- Constructor
    FluidFluid(const Foam::fvMesh& mesh);

    //- Configure
    bool configure(const YAML::Node adapterConfig);
    
    //- Add coupling data writers
    void addWriters(std::string dataName, Interface * interface);
    
    //- Add coupling data readers
    void addReaders(std::string dataName, Interface * interface);
};

}
}

#endif