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
	_error_("not implemented");
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
void HydrologyCuasAnalysis::UpdateWaterColumn(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		this->UpdateWaterColumn(element);
	}

}/*}}}*/
void HydrologyCuasAnalysis::UpdateWaterColumn(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating()) return;

	/*Intermediaries */
	IssmDouble  dt,drainage_rate,water_column;

	/*Retrieve all inputs and parameters*/
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);

	/*Get water column and drainage rate*/
	const int  numvertices= element->GetNumberOfVertices();
	IssmDouble* watercolumn  = xNew<IssmDouble>(numvertices);
	IssmDouble* drainagerate = xNew<IssmDouble>(numvertices);
	IssmDouble* meltingrate  = xNew<IssmDouble>(numvertices);
 	IssmDouble* watercolumn_max  = xNew<IssmDouble>(numvertices);
	element->GetInputListOnVertices(&watercolumn[0],WaterColumnOldEnum);
	element->GetInputListOnVertices(&drainagerate[0],HydrologyDrainageRateEnum);
	element->GetInputListOnVertices(&meltingrate[0],BasalforcingsGroundediceMeltingRateEnum);
	element->GetInputListOnVertices(&watercolumn_max[0],HydrologyWatercolumnMaxEnum);

	/*Add water*/
	for(int i=0;i<numvertices;i++){
		watercolumn[i] += (meltingrate[i]/rho_ice*rho_water-drainagerate[i])*dt;
		/*Check that water column height is within 0 and upper bound, correct if needed*/
 		if(watercolumn[i]>watercolumn_max[i]) watercolumn[i]=watercolumn_max[i];
 		if(watercolumn[i]<0) watercolumn[i]=0.;
	}

	/* Divide by connectivity, add degree of channelization as an input */
	/*FIXME: should be changed to P1, this is due to the NR, IsAllFloating will return 0 on this element, but it should not be DG*/
	element->AddInput(WatercolumnEnum,&watercolumn[0],P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(watercolumn);
	xDelete<IssmDouble>(meltingrate);
	xDelete<IssmDouble>(watercolumn_max);
	xDelete<IssmDouble>(drainagerate);
}/*}}}*/

void HydrologyCuasAnalysis::EffectiveTransmissivity(Element* element)/*{{{*/
	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating()) return;

	Input *S_input = element->GetInput(HydrologyStorageEnum); _assert_(S_input);
}/*}}}*/
void HydrologyCuasAnalysis::UpdateTransmissivity(Element* element){/*{{{*/
	/*
	 Update transmissivity...
	 */
	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating()) return;

	/*Intermediaries */
	IssmDouble  dt;

	/*Retrieve all inputs and parameters*/
	element->FindParam(&dt,HydrologySubTimeStepEnum);
	IssmDouble  rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);

	IssmDouble  latent    = element->FindParam(MaterialsLatentheatEnum);

	/*Get water column and drainage rate*/
	const int  numvertices= element->GetNumberOfVertices();
	IssmDouble* head         = xNew<IssmDouble>(numvertices);
	IssmDouble* T     = xNew<IssmDouble>(numvertices);
	IssmDouble* melt  = xNew<IssmDouble>(numvertices);
	IssmDouble* creep = xNew<IssmDouble>(numvertices);
	IssmDouble* cavity= xNew<IssmDouble>(numvertices);

	element->GetInputListOnVertices(&T[0],HydrologyCuasTrasmissivityEnum);
	element->GetInputListOnVertices(&creep[0],);
	element->GetInputListOnVertices(&melt[0],HydrologyCuasTrasmissivityEnum);
	element->GetInputListOnVertices(&cavity[0],HydrologyCuasTrasmissivityEnum);

	IssmDouble vel;
	Input *head_input = element->GetInput(HydrologyHeadEnum);      _assert_(head_input);
	Input *B_input    = element->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	Input *


	/* Divide by connectivity, add degree of channelization as an input */
	/*FIXME: should be changed to P1, this is due to the NR, IsAllFloating will return 0 on this element, but it should not be DG*/
	element->AddInput(WatercolumnEnum,&watercolumn[1],P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(watercolumn);
	xDelete<IssmDouble>(meltingrate);
	xDelete<IssmDouble>(watercolumn_max);
	xDelete<IssmDouble>(drainagerate);
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

