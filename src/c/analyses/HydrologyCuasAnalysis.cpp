#include "./HydrologyCuasAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyCuasAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	return;

}/*}}}*/
void HydrologyCuasAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	return;

}/*}}}*/
void HydrologyCuasAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	return;
}/*}}}*/
int  HydrologyCuasAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 0;
}/*}}}*/
void HydrologyCuasAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model,frictionlaw;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Cuas?*/
	if(hydrology_model!=HydrologycuasEnum) return;

	/*Add input to elements*/
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.conductivity",HydrologyChannelConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.Tinit",HydrologyCuasTinitEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.ss",HydrologyCuasSpecificStorageEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sy",HydrologyCuasSpecificYieldEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.roughness",HydrologyCuasRoughnessEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vel",VelEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);

	/*Update elements: */
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.Tinit",HydrologyCuasTrasmissivityEnum);
}/*}}}*/
void HydrologyCuasAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Cuas?*/
	if(hydrology_model!=HydrologycuasEnum) return;
	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ischannel_creep",HydrologyIschannelCreepEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ischannel_melt",HydrologyIschannelMeltEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ischannel_cavitiy",HydrologyIschannelCavityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.isconfined",HydrologyIsconfinedEnum));

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
	parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");

	/*Nothing else to add for now*/
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyCuasAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyCuasAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyCuasAnalysis::CreateDVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementMatrix* HydrologyCuasAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyCuasAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();
}/*}}}*/
ElementVector* HydrologyCuasAnalysis::CreatePVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyCuasAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyCuasAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyCuasAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyCuasAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/

/*Additional methods*/
void HydrologyCuasAnalysis::UpdateTransmissivity(Element* element){/*{{{*/
	/*
	 Update transmissivity...
	 */

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble* xyz_list=NULL;
	IssmDouble  dt;

	IssmDouble trans, trans_new;
	IssmDouble Pe; // effective pressure
	IssmDouble head;
	IssmDouble dh[3]; // deritavtives.
	IssmDouble dhdx, dhdy;
	IssmDouble n; /* Glen's flow exponent */
	IssmDouble B; /* ice rheology B */
	IssmDouble A; /* ice flow parameter A */

	IssmDouble bump_height, bump_space;

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,HydrologySubTimeStepEnum);
	IssmDouble  rho_ice   = basalelement->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = basalelement->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble  latent    = basalelement->FindParam(MaterialsLatentheatEnum);

	const int  numvertices= basalelement->GetNumberOfVertices();
	IssmDouble* trans_new = xNew<IssmDouble>(numvertices);

	IssmDouble vel;
	Input *vel_input  = basalelement->GetInput(VelEnum); _assert_(vel_input);
	Input *head_input = basalelement->GetInput(HydrologyHeadEnum);      _assert_(head_input);
	Input *B_input    = basalelement->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	Input *n_input    = basalelement->GetInput(MaterialsRheologyNEnum); _assert_(N_input);

	Input *conductivity_input=basalelement->GetInput(HydrologyConductivityEnum); _assert_(conductivity_input);
	Input *trans_input = basalelement->GetInput(HydrologyTransmissivityEnum); _assert_(trans_old_input);
	Input *Pe_input    = basalelement->GetInput(EffectivePressureEnum); _assert_(N_input);

	Input *bump_height_input = basalelement->GetInpu(HydrologyBumpHeightEnum); _assert_(bump_height_input);
	Input *bump_space_input  = basalelement->GetInpu(HydrologyBumpSpacingEnum); _assert_(bump_space_input);

	/* Update tarnsmissivity */
	Gauss *gauss;
	for(int iv=0;i<numvertices;i++){
		gauss->GaussVertex(iv);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		trans_old_input->GetInputValue(&trans,gauss);
		vel_input->GetInputValue(&vel,gauss);
		Pe_input->GetInputValue(&Pe,gauss);

		head_input->GetInputDerivativeValue(&dh[0],&xyz_list[0][0],&gauss);
		dhdx = dh[0];
		dhdy = dh[1];

		/*Get ice A parameters*/
		B_input->GetInputValue(&B, gauss);		
		n_input->GetInputValue(&n, gauss);
		A=pow(B,-n);

		/* geometry information */
		bump_height_input->GetInputValue(&bump_height,gauss);
		bump_space_input->GetInputValue(&bump_space,gauss);
		IssmDouble beta = bump_height/bump_space;

		/* Opening by melting */
		IssmDouble term1 = g * rho_water * trans / rho_ice / latent * (pow(dhdx,2.0) + pow(dhdy,2.0));
		IssmDouble term2 = 2*A*pow(n,-n)*pow(abs(N),n-1)*trans;
		IssmDouble term3 = beta*abs(vel)*conductivity;
		
		/* Resolve term in explictly */
		trans_new[iv] = dt*(term1 + term2 + term3) + trans 
	}

	/* Assign value */
	element->AddBasalInput(HydrologyTransmissivityNewEnum,trans_new,P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(trans_new);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void HyoldrogyCuasAnalysis::ComputeCavityOpening(Element* element){ /* {{{ */
	/*
	 Update cavity transmissivity...
	 */
	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating()) return;

	bool ischannel_cavitiy, ischannel_melt, ischannel_creep;
	IssmDouble dt;

	element->FindParam(&ischannel_cavitiy,HydrologyIschannelCavityEnum);
	element->FindParam(&ischannel_melt,HydrologyIschannelMeltEnum);
	element->FindParam(&ischannel_creep,HydrologyIschannelCreepEnum);
} /* }}} */
void HydrologyCuasAnalysis::UpdateStorage(Element* element){ /*{{{ */
	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	base_input=basalelement->GetInput(BaseEnum); _assert_(base_input);

	element->AddBasalInput(HydrologyStorageEnum,storage,P1DGEnum);
} /* }}} */
void HydroologyCuasAnalysis::GetEffectiveTransmissivity(Element *element){ /* {{{ */

} /* }}} */
