#ifndef FF_REDIAG_H
#define FF_REDIAG_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads velocity
class ReDiag : public CouplingDataUser
{

private:

    //- Re stresses field
    Foam::volSymmTensorField * ReStress_ ;
    
    //- Associated mesh
    const Foam::fvMesh& mesh_;
    
    //- Name of the Re stress field
    const std::string nameR_;

    //- Inlet volumetric flow rate (needed for mass correction)
    //double vDot_;
    
    //- Calculate mass correction for the values a the interface.
    //double massCorrection(double * buffer, int patchID);

public:

    //- Constructor
    ReDiag
    (
        const Foam::fvMesh& mesh,
        const std::string nameR
    );

    //- Write the velocity values into the buffer
    void write(double * buffer);

    //- Read the velocity values from the buffer
    void read(double * buffer);
};

}
}

// NOTE: In order to couple another variable, you need to add another
// subclass of the CouplingDataUser. Take this file and Velocity.C
// as an example and look for more notes in the Adapter.C and Adapter.H.
// Essentially, the new class needs to implement a constructor (which also
// sets the dataType_) and implement the write() and read() methods.
// Then, you need to create objects of this class in the Adapter.C,
// whenever requested from the configuration file (see notes).

#endif