#include "Velocity.H"
#include "cellSet.H"

using namespace Foam;

preciceAdapter::FF::Velocity::Velocity
(
    const Foam::fvMesh& mesh,
    const std::string nameU,
    const double vDot
)
:
U_(
    const_cast<volVectorField*>
    (
        &mesh.lookupObject<volVectorField>(nameU)
    )
),mesh_(mesh),vDot_(vDot)
{
    dataType_ = vector;
}

void preciceAdapter::FF::Velocity::write(double * buffer)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(U_->boundaryFieldRef()[patchID], i)
        {
            // Copy the velocity into the buffer
            // x-dimension
            buffer[bufferIndex++]
            =
            U_->boundaryFieldRef()[patchID][i].x();

            // y-dimension
            buffer[bufferIndex++]
            =
            U_->boundaryFieldRef()[patchID][i].y();

            // z-dimension
            buffer[bufferIndex++]
            =
            U_->boundaryFieldRef()[patchID][i].z();
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
			    buffer[bufferIndex++] = U_->internalField()[cells[i]].x();
			    buffer[bufferIndex++] = U_->internalField()[cells[i]].y();
			    buffer[bufferIndex++] = U_->internalField()[cells[i]].z();
        }
    }
}

double preciceAdapter::FF::Velocity::massCorrection(double * buffer, int patchID)
{
	// Calculate volumetric flow rate from the received data.
	double flowRate = 0;

	const polyPatch& cPatch = mesh_.boundaryMesh()[patchID];

	// For every cell of the patch
	//double area = 0;
	forAll(cPatch, i)
	{
		flowRate += buffer[3*i]*mesh_.magSf().boundaryField()[patchID][i];

		//area += mesh_->magSf().boundaryField()[patchID][i];
	}
 
  // Gather and add up the flowRate from all the processors
  reduce(flowRate,sumOp<double>());
  
  
	//std::cout << " Flow rate: " << flowRate  << " Area: " << area << std::endl;
	std::cout << " Flow rate: " << flowRate  << " Target flow rate: " << vDot_ << std::endl;

	return vDot_/flowRate;
}

void preciceAdapter::FF::Velocity::read(double * buffer)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);
        
        // Mass correction
        // TODO: I'm not sure how this will work with more than one patch, specially patches that are not in the x direction or don't cover the whole yz plane. 
        // Also, this only works if the flow direction perpendicular to the patch is ux. 
		    double mcorrection = 1.0;
		    if(vDot_ > 0)
			    mcorrection = massCorrection(buffer, patchID);
                          
        std::cout << " Mass correction: " << mcorrection << std::endl;

        // For every cell of the patch
        forAll(U_->boundaryFieldRef()[patchID], i)
        {
            // Set the velocity as the buffer value
            // x-dimension
            U_->boundaryFieldRef()[patchID][i].x()
            =
            mcorrection * buffer[bufferIndex++];

            // y-dimension
            U_->boundaryFieldRef()[patchID][i].y()
            =
            buffer[bufferIndex++];

            // z-dimension
            U_->boundaryFieldRef()[patchID][i].z()
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
			    U_->ref()[cells[i]].x() = buffer[bufferIndex++];
			    U_->ref()[cells[i]].y() = buffer[bufferIndex++];
			    U_->ref()[cells[i]].z() = buffer[bufferIndex++];
        }
    }
    
    
}