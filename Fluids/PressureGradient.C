#include "PressureGradient.H"
#include "cellSet.H"
#include "fixedGradientFvPatchField.H"

using namespace Foam;

preciceAdapter::Fluids::PressureGradient::PressureGradient
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

void preciceAdapter::Fluids::PressureGradient::write(double * buffer)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++ )
    {
        int patchID = patchIDs_.at(j);
        
        P_->boundaryFieldRef()[patchID].updateCoeffs();
        
        fixedGradientFvPatchField<double>* fixedGradP = dynamic_cast<fixedGradientFvPatchField<double>*>(&P_->boundaryFieldRef()[patchID]);
        if (fixedGradP == nullptr)
        {
          std::cout << "Dynamic cast in write method for pressureGradient failed. Please make sure that the pressure BC in the interface is set to fixedGradient." << std::endl; 
        }
        
        // For every cell of the patch
        for( int i=0; i < fixedGradP->gradient().size(); i++)
        {
            // Copy the pressure gradient into the buffer
             buffer[bufferIndex++] 
             =
             fixedGradP->gradient()[i];            
        }
                
    }
    
    for (uint j = 0 ; j < cellSetNames_.size(); j++)
    {
        // TODO: Caclulate the pressure gradient at the cells that form the cellSet and write it to the buffer. 
        // However, which direction should the gradient be in? cellSets doesn't have normal vector. 
        std::cout << " Writing pressure gradient from cellSets to buffer is not yet implemented. Pressure gradient information is not copied to the buffer." << std::endl;
    }
}

void preciceAdapter::Fluids::PressureGradient::read(double * buffer)
{

	int bufferIndex = 0;

 	std::cout << "Reading pressure from LUMIS: " << std::endl;

	// For every boundary patch of the interface
	for (uint j = 0; j < patchIDs_.size(); j++)
	{
		int patchID = patchIDs_.at(j);

		P_->boundaryFieldRef()[patchID].updateCoeffs();
		Field<double> grad = P_->boundaryFieldRef()[patchID];
        // For every cell of the patch
    //forAll(grad, i)
    for( int i=0; i < grad.size(); i++)
    {
      // Copy the pressure gradient into grad
         grad[i] = buffer[bufferIndex++];            
    }
        
    fixedGradientFvPatchField<double>* fixedGradP = dynamic_cast<fixedGradientFvPatchField<double>*>(&P_->boundaryFieldRef()[patchID]);
    if (fixedGradP == nullptr)
    {
      std::cout << "Dynamic cast in read method for pressure failed. Please make sure that the pressure BC in the interface is set to fixedGradient." << std::endl; 
    }
    //std::cout << "Gradient field extracted: " << std::endl;
   
    fixedGradP->gradient() = grad; 
    //std::cout << "gradient field set: " << std::endl;

	}
 
  for (uint j = 0 ; j < cellSetNames_.size(); j++)
  {
      // TODO: Caclulate the pressure gradient at the cells that form the cellSet and write it to the buffer. 
      // However, which direction should the gradient be in? cellSets doesn't have normal vector. 
      std::cout << " Writing pressure gradient from cellSets to buffer is not yet implemented. Pressure gradient information is not copied to the buffer." << std::endl;
  }

}
