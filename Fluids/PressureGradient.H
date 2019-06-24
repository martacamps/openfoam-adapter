#ifndef FLUIDS_PRESSUREGRADIENT_H
#define FLUIDS_PRESSUREGRADIENT_H

#include "CouplingDataUser.H"
#include "fvMesh.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace Fluids
{

//- Class that writes and reads pressure.
class PressureGradient : public CouplingDataUser
{

private:

    //- Pressure field
    Foam::volScalarField * P_;

    //- Associated mesh
    const Foam::fvMesh& mesh_;

public:

    //- Constructor
    PressureGradient
    (
        const Foam::fvMesh& mesh,
        const std::string nameT
    );

    //- Write the pressure values into the buffer
    void write(double * buffer);

    //- Read the pressure values from the buffer
    void read(double * buffer);

};

}
}

// NOTE: In order to couple another variable, you need to add another
// subclass of the CouplingDataUser. Take this file and PressureGradient.C
// as an example and look for more notes in the Adapter.C and Adapter.H.
// Essentially, the new class needs to implement a constructor (which also
// sets the dataType_) and implement the write() and read() methods.
// Then, you need to create objects of this class in the Adapter.C,
// whenever requested from the configuration file (see notes).

#endif
