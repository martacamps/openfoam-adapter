#include "Epsilon.H"
#include "cellSet.H"
// Include some type of turbulence model .H?

using namespace Foam;

preciceAdapter::FF::Epsilon::Epsilon
(
    const Foam::fvMesh& mesh,
    const std::string nameE
)
:
epsilon_(NULL),mesh_(mesh),nameE_(nameE)
{
    dataType_ = scalar; 
}

void preciceAdapter::FF::Epsilon::write(double * buffer)
{
    int bufferIndex = 0;
    
    // Get a pointer to the epsilon field if it didn't exist
    if (epsilon_ == NULL)
      epsilon_ = const_cast<volScalarField*>(&mesh_.lookupObject<volScalarField>(nameE_));

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(epsilon_->boundaryFieldRef()[patchID], i)
        {
            // Copy the pressure into the buffer
            buffer[bufferIndex++]
            =
            epsilon_->boundaryFieldRef()[patchID][i];
        }
    }
    
    // For every cell set of the interface
    // TODO: Do I have to create the cellSet each time? Don't they have indices or something like patches do?
		// Maybe I can store pointers to the cellSets?
    for (uint j = 0; j < cellSetNames_.size(); j++)
    {
        cellSet overlapRegion(mesh_, cellSetNames_.at(j));
		    const labelList & cells = overlapRegion.toc();
        for( uint i=0; i < cells.size(); i++)
        {
        	// Copy the pressure into the buffer
        	buffer[bufferIndex++]
        	=
        	epsilon_->internalField()[cells[i]];
        }
    }
}

void preciceAdapter::FF::Epsilon::read(double * buffer)
{
    int bufferIndex = 0;
    
    // Get a pointer to the epsilon field if it didn't exist
    if (epsilon_ == NULL)
      epsilon_ = const_cast<volScalarField*>(&mesh_.lookupObject<volScalarField>(nameE_));

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(epsilon_->boundaryFieldRef()[patchID], i)
        {
            // Set the pressure as the buffer value
            epsilon_->boundaryFieldRef()[patchID][i]
            =
            buffer[bufferIndex++];
        }
    }
    
     // For every cell set of the interface
    // TODO: Do I have to create the cellSet each time? Don't they have indices or something like patches do?
		// Maybe I can store pointers to the cellSets?
    for (uint j = 0; j < cellSetNames_.size(); j++)
    {
        cellSet overlapRegion(mesh_, cellSetNames_.at(j));
		    const labelList & cells = overlapRegion.toc();
        for( uint i=0; i < cells.size(); i++)
        {
        	// Set the pressure as the buffer value
          epsilon_->ref()[cells[i]]
          = 
          buffer[bufferIndex++];
        }
    }
}