#ifndef FF_EPSILON_H
#define FF_EPSILON_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads pressure
class Epsilon : public CouplingDataUser
{

private:

    //- Epsilon field
    Foam::volScalarField * epsilon_;
    
    //- Associated mesh
    const Foam::fvMesh& mesh_;

public:

    //- Constructor
    Epsilon
    (
        const Foam::fvMesh& mesh,
        const std::string nameE
    );

    //- Write the pressure values into the buffer
    void write(double * buffer);

    //- Read the pressure values from the buffer
    void read(double * buffer);
};

}
}

// NOTE: In order to couple another variable, you need to add another
// subclass of the CouplingDataUser. Take this file and Epsilon.C
// as an example and look for more notes in the Adapter.C and Adapter.H.
// Essentially, the new class needs to implement a constructor (which also
// sets the dataType_) and implement the write() and read() methods.
// Then, you need to create objects of this class in the Adapter.C,
// whenever requested from the configuration file (see notes).

#endif