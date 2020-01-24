#include "Pressure.H"
#include "cellSet.H"

using namespace Foam;

preciceAdapter::FF::Pressure::Pressure
(
    const Foam::fvMesh& mesh,
    const std::string nameP
)
:
p_(
    const_cast<volScalarField*>
    (
        &mesh.lookupObject<volScalarField>(nameP)
    )
),mesh_(mesh)
{
    dataType_ = scalar;
}

void preciceAdapter::FF::Pressure::write(double * buffer)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(p_->boundaryFieldRef()[patchID], i)
        {
            // Copy the pressure into the buffer
            buffer[bufferIndex++]
            =
            p_->boundaryFieldRef()[patchID][i];
        }
    }
    
    // For every cell set of the interface
    // TODO: Do I have to create the cellSet each time? Don't they have indices or something like patches do?
		// Maybe I can store pointers to the cellSets?
    for (uint j = 0; j < cellSetNames_.at(j); j++)
    {
        cellSet overlapRegion(mesh_, cellSetNames_.at(j));
		    const labelList & cells = overlapRegion.toc();
        for( uint i=0; i < cells.size(); i++)
        {
        	// Copy the pressure into the buffer
        	buffer[bufferIndex++]
        	=
        	P_->internalField()[cells[i]];
        }
    }
}

void preciceAdapter::FF::Pressure::read(double * buffer)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(p_->boundaryFieldRef()[patchID], i)
        {
            // Set the pressure as the buffer value
            p_->boundaryFieldRef()[patchID][i]
            =
            buffer[bufferIndex++];
        }
    }
    
     // For every cell set of the interface
    // TODO: Do I have to create the cellSet each time? Don't they have indices or something like patches do?
		// Maybe I can store pointers to the cellSets?
    for (uint j = 0; j < cellSetNames_.at(j); j++)
    {
        cellSet overlapRegion(mesh_, cellSetNames_.at(j));
		    const labelList & cells = overlapRegion.toc();
        for( uint i=0; i < cells.size(); i++)
        {
        	// Set the pressure as the buffer value
          p_->ref()[cells[i]].x()
          = 
          buffer[bufferIndex++];
        }
    }
}