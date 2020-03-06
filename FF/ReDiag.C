#include "ReDiag.H"
#include "cellSet.H"

using namespace Foam;

preciceAdapter::FF::ReDiag::ReDiag
(
    const Foam::fvMesh& mesh,
    const std::string nameR
)
:
/*ReStress_(
    const_cast<volSymmTensorField*>
    (
        &mesh.lookupObject<volSymmTensorField>(nameR)
    )
)*/nameR_(nameR),mesh_(mesh),ReStress_(NULL)
{
    dataType_ = vector;
}

void preciceAdapter::FF::ReDiag::write(double * buffer)
{
    int bufferIndex = 0;
    
    // Get the Reynolds stress 
    //const incompressible::turbulenceModel& model = &mesh_.lookupObject<incompressible::turbulenceModel>(turbulenceModel::propertiesName);
    
    //volSymmTensorField& ReStress = model.R();
    
    // Get a pointer to the Re stress if it didn't exist
    if (ReStress_ == NULL)
      ReStress_ = const_cast<volSymmTensorField*>(&mesh_.lookupObject<volSymmTensorField>(nameR_));

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(ReStress_->boundaryFieldRef()[patchID], i)
        {
            // Copy the velocity into the buffer
            // x-dimension
            buffer[bufferIndex++]
            =
            ReStress_->boundaryFieldRef()[patchID][i].xx();

            // y-dimension
            buffer[bufferIndex++]
            =
            ReStress_->boundaryFieldRef()[patchID][i].yy();

            // z-dimension
            buffer[bufferIndex++]
            =
            ReStress_->boundaryFieldRef()[patchID][i].zz();
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
        	// Copy the three components of the velocity into the buffer
			    buffer[bufferIndex++] = ReStress_->internalField()[cells[i]].xx();
			    buffer[bufferIndex++] = ReStress_->internalField()[cells[i]].yy();
			    buffer[bufferIndex++] = ReStress_->internalField()[cells[i]].zz();
        }
    }
}

void preciceAdapter::FF::ReDiag::read(double * buffer)
{
    int bufferIndex = 0;
    
    // Get the Reynolds stress 
    //const incompressible::turbulenceModel& model = &mesh_.lookupObject<incompressible::turbulenceModel>(turbulenceModel::propertiesName);
    //volSymmTensorField& ReStress = model.R();
    
    // Get a pointer to the Re stress if it didn't exist
    if (ReStress_ == NULL)
      ReStress_ = const_cast<volSymmTensorField*>(&mesh_.lookupObject<volSymmTensorField>(nameR_));

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);
        
                         
        // For every cell of the patch
        forAll(ReStress_->boundaryFieldRef()[patchID], i)
        {
            // Set the velocity as the buffer value
            // x-dimension
            ReStress_->boundaryFieldRef()[patchID][i].xx()
            =
            buffer[bufferIndex++];

            // y-dimension
            ReStress_->boundaryFieldRef()[patchID][i].yy()
            =
            buffer[bufferIndex++];

            // z-dimension
            ReStress_->boundaryFieldRef()[patchID][i].zz()
            =
            buffer[bufferIndex++];
        }
    }
    
    // For every cell set of the interface
    // TODO: Do I have to create the cellSet each time? Don't they have indices or something like patches do?
		// Maybe I can store pointers to the cellSets?
    // TODO: Implement mass correction for cellSets?
    for (uint j = 0; j < cellSetNames_.size(); j++)
    {
        cellSet overlapRegion(mesh_, cellSetNames_.at(j));
		    const labelList & cells = overlapRegion.toc();
        for( uint i=0; i < cells.size(); i++)
        {
        	// Set the velocity to the buffer value
			    ReStress_->ref()[cells[i]].xx() = buffer[bufferIndex++];
			    ReStress_->ref()[cells[i]].yy() = buffer[bufferIndex++];
			    ReStress_->ref()[cells[i]].zz() = buffer[bufferIndex++];
        }
    }
    
    
}