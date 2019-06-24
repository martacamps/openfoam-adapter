#include "Pressure.H"
#include "cellSet.H"
//#include "readGradientFvPatchField.H"
#include "fixedGradientFvPatchField.H"

using namespace Foam;

preciceAdapter::Fluids::Pressure::Pressure
(
    const Foam::fvMesh& mesh,
    const std::string nameT
)
:
P_(
    const_cast<volScalarField*>
    (
        &mesh.lookupObject<volScalarField>(nameT)
    )
),mesh_(mesh)
{
    dataType_ = scalar;
}

void preciceAdapter::Fluids::Pressure::write(double * buffer)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++ )
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(P_->boundaryFieldRef()[patchID], i)
        {
            // Copy the pressure into the buffer
            buffer[bufferIndex++]
            =
            P_->boundaryFieldRef()[patchID][i];
        }
        
    }

    // For every cellSet of the interface
    for (uint j = 0; j < cellSetNames_.size(); j++)
    {
        // TODO: Do I have to create the cellSet each time? Don't they have indices or something like patches do?
	    	// Maybe I can store pointers to the cellSets?
		    cellSet overlapRegion(mesh_, cellSetNames_[j]);
		    const labelList & cells = overlapRegion.toc();
        for( int i=0; i < cells.size(); i++)
        {
        	// Copy the pressure into the buffer
        	buffer[bufferIndex++]
        	=
        	P_->internalField()[cells[i]];
        }
    }
}

void preciceAdapter::Fluids::Pressure::read(double * buffer)
{

	int bufferIndex = 0;

 	std::cout << "Reading pressure from LUMIS: " << std::endl;

	// For every boundary patch of the interface
	for (uint j = 0; j < patchIDs_.size(); j++)
	{
		int patchID = patchIDs_.at(j);
   
    // For every cell of the patch
    forAll(P_->boundaryFieldRef()[patchID], i)
    {
      // Set the pressure as the buffer value
      P_->boundaryFieldRef()[patchID][i]
      =
      buffer[bufferIndex++];
    }
	}
 
  // For every cellSet of the interface
  for (uint j = 0; j < cellSetNames_.size(); j++)
  {
      // TODO: Do the cellSets have to be created before using all the time ? Maybe I can store pointers to the cellSets?
	    cellSet overlapRegion(mesh_, cellSetNames_.at(j));
	    const labelList & cells = overlapRegion.toc();
             
  		for( int i=0; i < cells.size(); i++)
  		{
  			// Set the velocity to the buffer value
  			P_->ref()[cells[i]] = buffer[bufferIndex++];
  		}
  }

}
