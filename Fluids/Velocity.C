#include "Velocity.H"
#include "cellSet.H"
//#include "vector.H"

using namespace Foam;

preciceAdapter::Fluids::Velocity::Velocity
(
    const Foam::fvMesh& mesh,
    const std::string nameT,
    const double vDot
)
:
U_(
    const_cast<volVectorField*>
    (
        &mesh.lookupObject<volVectorField>(nameT)
    )
),mesh_(mesh),vDot_(vDot)
{
    dataType_ = vector;
}

void preciceAdapter::Fluids::Velocity::write(double * buffer)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++ )
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(U_->boundaryFieldRef()[patchID], i)
        {
            // Copy the three components of the velocity into the buffer
            buffer[bufferIndex++] = U_->boundaryFieldRef()[patchID][i].x();
            buffer[bufferIndex++] = U_->boundaryFieldRef()[patchID][i].y();
            buffer[bufferIndex++] = U_->boundaryFieldRef()[patchID][i].z();
        }
        
    }
    
    for (uint j = 0; j < cellSetNames_.size(); j++)
    {
        // TODO: Do I have to create the cellSet each time? Don't they have indices or something like patches do?
        // Maybe I can store pointers to the cellSets?
        cellSet overlapRegion(mesh_, cellSetNames_.at(j));
		    const labelList & cells = overlapRegion.toc();
		    for( int i=0; i < cells.size(); i++)
		    {
			      // Copy the three components of the velocity into the buffer
			      buffer[bufferIndex++] = U_->internalField()[cells[i]].x();
			      buffer[bufferIndex++] = U_->internalField()[cells[i]].y();
			      buffer[bufferIndex++] = U_->internalField()[cells[i]].z();
		    }
    }

}

double preciceAdapter::Fluids::Velocity::massCorrection(double * buffer, int patchID)
{

	// Calculate volumetric flow rate from the LUMIS data.
	double flowRate = 0;

	const polyPatch& cPatch = mesh_.boundaryMesh()[patchID];

	// For every cell of the patch
	//double area = 0;
  // TODO: Use the patch normal vector so that massCorrection is valid for other than patches perpendicular to x direction. 
	forAll(cPatch, i)
	{
     //Foam::vector velocity(buffer[3*i], buffer[3*i+1], buffer[3*i+2]);
     
     //flowRate += (velocity & cPatch.faceAreas()[i]) * mesh_.magSf().boundaryField()[patchID][i];
 
		 flowRate += buffer[3*i]*mesh_.magSf().boundaryField()[patchID][i];

		//area += mesh_->magSf().boundaryField()[patchID][i];
	}
 
  // Gather and add up the flowRate from all the processors
  reduce(flowRate,sumOp<double>());
  
  
	//std::cout << " Flow rate: " << flowRate  << " Area: " << area << std::endl;
	std::cout << " Flow rate: " << flowRate  << " Target flow rate: " << vDot_ << std::endl;

	return vDot_/flowRate;
}

void preciceAdapter::Fluids::Velocity::read(double * buffer)
{

    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Mass correction
        // It assumes that all the patches have the same flow rate vDot_
        // Also that all the patches go completelly across the domain. 
		    float mcorrection = 1.0;
		    if(vDot_ > 0)
        {
		        mcorrection = massCorrection(buffer, patchID);
        }

		    std::cout << " Mass correction: " << mcorrection << std::endl;

        // For every cell of the patch
        // TODO: n has to be the patch normal instead of just 1 0 0. 
        Foam::vector n(1.0, 0.0, 0.0);
        forAll(U_->boundaryFieldRef()[patchID], i)
        {
            // Get normal vector 
            //Foam::vector n = mesh_.boundaryMesh()[patchID].faceAreas()[i];
        
            // Set the velocity as the buffer value
            U_->boundaryFieldRef()[patchID][i].x() = mcorrection * n.x() * buffer[bufferIndex++];
            U_->boundaryFieldRef()[patchID][i].y() = mcorrection * n.y() * buffer[bufferIndex++];
            U_->boundaryFieldRef()[patchID][i].z() = mcorrection * n.z() * buffer[bufferIndex++];

            // Check that the mass flow has been corrected?
        }
    }

    //TODO: Cacluate and apply mass correction to cellSets
    for (uint j = 0; j < cellSetNames_.size(); j++)
    {
        // TODO: Do the cellSets have to be created before using all the time ? Maybe I can store pointers to the cellSets?
		    cellSet overlapRegion(mesh_, cellSetNames_.at(j));
		    const labelList & cells = overlapRegion.toc();
            
        if(vDot_ > 0)
        {
		        std::cout << " Mass correction not implemented for cellSets." << std::endl;
        }
        
    		for( int i=0; i < cells.size(); i++)
    		{
    			// Set the velocity to the buffer value
    			U_->ref()[cells[i]].x() = buffer[bufferIndex++];
    			U_->ref()[cells[i]].y() = buffer[bufferIndex++];
    			U_->ref()[cells[i]].z() = buffer[bufferIndex++];
    		}
    }
}
