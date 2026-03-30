#include "./HydrologyCuasAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

#define FINITEELEMENT P1Enum

/*Model processing*/
void HydrologyCuasAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologycuasEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spchead",HydrologyCuasAnalysisEnum,P1Enum);
	
}/*}}}*/
void HydrologyCuasAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Cuas?*/
	if(hydrology_model!=HydrologycuasEnum) return;

	//NOTE: Reuse modules in HydrologyCuasAnalysis
	/*Create discrete loads for Moulins*/
	_printf0_("HydrologyCuasAnalysis::CreateLoads: Add moulins for inputs\n");
	CreateSingleNodeToElementConnectivity(iomodel);
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(1,"md.mesh.vertexonbase");
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}	
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");

	/*Deal with Neumann BC*/
	int M,N;
	int *segments = NULL;
	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(&segments,&M,&N,"md.mesh.segments2d");
	}
	else if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchData(&segments,&M,&N,"md.mesh.segments");
	}
	else{
		_error_("mesh type not supported yet");
	}

	/*Check that the size seem right*/
	_assert_(N==3); _assert_(M>=3);

	for(int i=0;i<M;i++){
		if(iomodel->my_elements[segments[i*3+2]-1]){
			loads->AddObject(new Neumannflux(i+1,i,iomodel,segments));
		}
	}
	xDelete<int>(segments);

}/*}}}*/
void HydrologyCuasAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want CUAS?*/
	if(hydrology_model!=HydrologycuasEnum) return;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,HydrologyCuasAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");

}/*}}}*/
int  HydrologyCuasAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyCuasAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model,frictionlaw;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Cuas?*/
	if(hydrology_model!=HydrologycuasEnum) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Add input to elements*/
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.geothermalflux",BasalforcingsGeothermalfluxEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxBaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyBaseEnum);
	}
	
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.head",HydrologyHeadEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.conductivity",HydrologyConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.layer_thickness",HydrologySheetThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.transmissivity",HydrologyTransmissivityEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.ss",HydrologySpecificStorageEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sy",HydrologySpecificYieldEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.englacial_input",HydrologyEnglacialInputEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_spacing",HydrologyBumpSpacingEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.moulin_input",HydrologyMoulinInputEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.neumannflux",HydrologyNeumannfluxEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);


	/*Initialize requested outputs in case they are not defined later for this partition*/
	iomodel->ConstantToInput(inputs,elements,0.,HydrologyChannelMeltRateEnum,FINITEELEMENT);
	iomodel->ConstantToInput(inputs,elements,0.,HydrologyChannelCreepRateEnum,FINITEELEMENT);
	iomodel->ConstantToInput(inputs,elements,0.,HydrologyChannelCavityRateEnum,FINITEELEMENT);
	iomodel->ConstantToInput(inputs,elements,0.,HydrologyBasalFluxEnum,FINITEELEMENT);
	iomodel->ConstantToInput(inputs,elements,0.,HydrologyWaterVxEnum,FINITEELEMENT);
	iomodel->ConstantToInput(inputs,elements,0.,HydrologyWaterVyEnum,FINITEELEMENT);
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
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ischannel_cavity",HydrologyIschannelCavityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.isconfined",HydrologyIsConfinedEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.melt_flag",HydrologyMeltFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.unconfinedSmooth",HydrologyCuasUnconfinedSmoothEnum));

	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.Tmax",HydrologyCuasTmaxEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.Tmin",HydrologyCuasTminEnum));

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
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyCuasAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyCuasAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;

	IssmDouble storage, Teff;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke  = basalelement->NewElementMatrix();
	IssmDouble* dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble* basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);

	/*Get Params*/
	IssmDouble dt;
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);

	/*Get all inputs and parameters*/
	Input* storage_input = basalelement->GetInput(HydrologyStorageEnum); _assert_(storage_input);
	Input* Teff_input	= basalelement->GetInput(HydrologyTransmissivityEffectiveEnum); _assert_(Teff_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(3);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		storage_input->GetInputValue(&storage,gauss);
		Teff_input->GetInputValue(&Teff,gauss);

		/* Diffusion term */
		IssmDouble factor1 = Teff*gauss->weight*Jdet;
		
		/* Transient mass term */
		IssmDouble factor2 = storage/dt*gauss->weight*Jdet;

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += factor1*(dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j])
				  + factor2*basis[i]*basis[j];
			}
		}
	}

	/* Clean up and return */
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* HydrologyCuasAnalysis::CreatePVector(Element* element){/*{{{*/
	
	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  headold,storage;
	IssmDouble  meltrate;
	IssmDouble  englacial_rate;
	IssmDouble  rho_water, rho_ice;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe	= basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	Input* headold_input	 = basalelement->GetInput(HydrologyHeadOldEnum); _assert_(headold_input);
	Input* storage_input  = basalelement->GetInput(HydrologyStorageEnum); _assert_(storage_input);
	Input* meltrate_input = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(meltrate_input);
	Input* englacial_input = basalelement->GetInput(HydrologyEnglacialInputEnum); _assert_(englacial_input);

	/*Get Params*/
	IssmDouble dt;
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&rho_water,MaterialsRhoFreshwaterEnum);
	basalelement->FindParam(&rho_ice,MaterialsRhoIceEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(3);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		headold_input->GetInputValue(&headold,gauss);
		storage_input->GetInputValue(&storage,gauss);
		meltrate_input->GetInputValue(&meltrate,gauss);
		englacial_input->GetInputValue(&englacial_rate,gauss);

		meltrate = meltrate*rho_ice/rho_water;
		IssmDouble factor = Jdet*gauss->weight*
		  ((meltrate + englacial_rate) + storage*headold/dt);

		for(int i=0;i<numnodes;i++){
			pe->values[i] += factor*basis[i];
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void           HydrologyCuasAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,HydrologyHeadEnum);
}/*}}}*/
void           HydrologyCuasAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyCuasAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	 /*Only update if on base*/
	if(!element->IsOnBase()) return;

	/*Intermediary*/
	IssmDouble dh[3];
	int* doflist = NULL;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];
		
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Add input to the element: */
	element->AddBasalInput(HydrologyHeadEnum,values,element->GetElementType());


	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);

}/*}}}*/
void           HydrologyCuasAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Update active elements based on ice levelset and ocean levelset*/
	GetMaskOfIceVerticesLSMx(femmodel,true);
	SetActiveNodesLSMx(femmodel,true);

	IssmDouble rho_ice   = femmodel->parameters->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = femmodel->parameters->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = femmodel->parameters->FindParam(ConstantsGEnum);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){

		/*Get current element and return if not on base*/
		Element *element  = xDynamicCast<Element*>(object);
		if(!element->IsOnBase()) continue;

		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *bed       = xNew<IssmDouble>(numnodes);
		IssmDouble *thickness = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&bed[0],BaseEnum);
		element->GetInputListOnNodes(&thickness[0],ThicknessEnum);
		element->GetInputListOnNodes(&ls_active[0],HydrologyMaskNodeActivationEnum);

		//for(int in=0;in<numnodes;in++){ //
		for(int in=0;in<3;in++){ //
			Node* node=element->GetNode(in);
			if(mask[in]>0. && ls_active[in]==1.){
				node->Activate(); //Not sure if we need this!
			}
			else{
				node->Deactivate();// node should be inactive
				node->ApplyConstraint(0,0.0); // set head (dof 0) to 0.0 m
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(bed);
		xDelete<IssmDouble>(thickness);
		xDelete<IssmDouble>(ls_active);
	}

	return;
}/*}}}*/

/*Additional methods*/
void HydrologyCuasAnalysis::UpdateChannelRates(FemModel* femmodel){/*{{{*/
	 for(Object* & object : femmodel->elements->objects){
		  Element *element = xDynamicCast<Element*>(object);
		  UpdateChannelRates(element);
	 }
}/*}}}*/
void HydrologyCuasAnalysis::UpdateChannelRates(Element* element){/*{{{*/
	/*
	 Update channel rates
	 1) melt opening
	 2) creep closure
	 3) cavity opening
	 */

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble* xyz_list=NULL;
	IssmDouble  dt;
	IssmDouble  Jdet;

	IssmDouble  vel,vx,vy;
	IssmDouble  trans;
	IssmDouble  conductivity;
	IssmDouble  Neff; // effective pressure
	IssmDouble  head;
	IssmDouble  dh[2]; // deritavtives.
	IssmDouble  dhdx, dhdy;
	IssmDouble  n; /* Glen's flow exponent */
	IssmDouble  B; /* ice rheology B */
	IssmDouble  A; /* ice flow parameter A */

	IssmDouble bump_height, bump_space;

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  g         = basalelement->FindParam(ConstantsGEnum);
	IssmDouble  rho_ice   = basalelement->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = basalelement->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble  latent    = basalelement->FindParam(MaterialsLatentheatEnum);
	bool ischannel_melt, ischannel_creep, ischannel_cavity;
	basalelement->FindParam(&ischannel_creep,HydrologyIschannelCreepEnum);
	basalelement->FindParam(&ischannel_melt,HydrologyIschannelMeltEnum);
	basalelement->FindParam(&ischannel_cavity,HydrologyIschannelCavityEnum);

	const int  numvertices= basalelement->GetNumberOfVertices();
	IssmDouble *channel_creep, *channel_melt, *channel_cavity;
	channel_melt  = xNew<IssmDouble>(numvertices);
	channel_creep = xNew<IssmDouble>(numvertices);
	channel_cavity= xNew<IssmDouble>(numvertices);

	Input *vx_input  = basalelement->GetInput(VxEnum); _assert_(vx_input);
	Input *vy_input  = basalelement->GetInput(VxEnum); _assert_(vy_input);
	Input *head_input = basalelement->GetInput(HydrologyHeadEnum);      _assert_(head_input);
	Input *B_input    = basalelement->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	Input *n_input    = basalelement->GetInput(MaterialsRheologyNEnum); _assert_(n_input);

	Input *conductivity_input=basalelement->GetInput(HydrologyConductivityEnum); _assert_(conductivity_input);
	Input *trans_input = basalelement->GetInput(HydrologyTransmissivityOldEnum); _assert_(trans_input);
	Input *Neff_input  = basalelement->GetInput(EffectivePressureEnum); _assert_(Neff_input);

	Input *bump_height_input = basalelement->GetInput(HydrologyBumpHeightEnum); _assert_(bump_height_input);
	Input *bump_space_input  = basalelement->GetInput(HydrologyBumpSpacingEnum); _assert_(bump_space_input);

	/* Update transmissivity */
	Gauss *gauss=basalelement->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussNode(basalelement->GetElementType(),iv);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		trans_input->GetInputValue(&trans,gauss);
		conductivity_input->GetInputValue(&conductivity,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel = sqrt(vx*vx + vy*vy);
		Neff_input->GetInputValue(&Neff,gauss);

		head_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
		dhdx = dh[0];
		dhdy = dh[1];

		/* Get ice A parameters*/
		B_input->GetInputValue(&B, gauss);		
		n_input->GetInputValue(&n, gauss);
		A=pow(B,-n);

		/* geometry information */
		bump_height_input->GetInputValue(&bump_height,gauss);
		bump_space_input->GetInputValue(&bump_space,gauss);
		IssmDouble beta = bump_height/bump_space;

		/* Opening by melting */
		if(ischannel_melt)   channel_melt[iv]= g*rho_water*conductivity*trans/rho_ice/latent*(pow(dhdx,2.0) + pow(dhdy,2.0));
		else channel_creep[iv] = 0.0;

		/* Closure by creep */
		if(ischannel_creep)  channel_creep[iv]= -2*A*pow(n,-n)*pow(abs(Neff),n-1)*Neff*trans;
		else channel_creep[iv] = 0.0;

		/* Opening by cavity */
		if(ischannel_cavity) channel_cavity[iv]= beta*abs(vel)*conductivity;
		else channel_cavity[iv] = 0.0;
	}

	/* Assign value */
	element->AddBasalInput(HydrologyChannelCavityRateEnum,channel_cavity,FINITEELEMENT);
	element->AddBasalInput(HydrologyChannelMeltRateEnum,channel_melt,FINITEELEMENT);
	element->AddBasalInput(HydrologyChannelCreepRateEnum,channel_creep,FINITEELEMENT);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(channel_melt);
	xDelete<IssmDouble>(channel_creep);
	xDelete<IssmDouble>(channel_cavity);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/

void HydrologyCuasAnalysis::UpdateTransmissivity(FemModel* femmodel){/*{{{*/
	 for(Object* & object : femmodel->elements->objects){
		  Element *element = xDynamicCast<Element*>(object);
		  UpdateTransmissivity(element);
	 }
}/*}}}*/
void HydrologyCuasAnalysis::UpdateTransmissivity(Element* element){/*{{{*/
	/*
	 Update transmissivity using updated creep, melt, cavity opening
	 */

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble* xyz_list=NULL;
	IssmDouble  dt;
	IssmDouble  Jdet;

	IssmDouble  trans;
	IssmDouble* trans_new;
	IssmDouble  channel_melt, channel_cavity, channel_creep;

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  Tmin = basalelement->FindParam(HydrologyCuasTminEnum);
	IssmDouble  Tmax = basalelement->FindParam(HydrologyCuasTmaxEnum);

	const int  numvertices= basalelement->GetNumberOfVertices();
	trans_new = xNew<IssmDouble>(numvertices);

	Input* trans_input       = basalelement->GetInput(HydrologyTransmissivityOldEnum); _assert_(trans_input);
	Input* channelmelt_input = basalelement->GetInput(HydrologyChannelMeltRateEnum); _assert_(channelmelt_input);
	Input* channelcavity_input = basalelement->GetInput(HydrologyChannelCavityRateEnum); _assert_(channelcavity_input);
	Input* channelcreep_input = basalelement->GetInput(HydrologyChannelCreepRateEnum); _assert_(channelcreep_input);

	/* Update transmissivity */
	Gauss *gauss=basalelement->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussNode(basalelement->GetElementType(),iv);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		trans_input->GetInputValue(&trans,gauss);

		channelmelt_input->GetInputValue(&channel_melt, gauss);
		channelcavity_input->GetInputValue(&channel_cavity, gauss);
		channelcreep_input->GetInputValue(&channel_creep, gauss);
		
		/* Resolve term in explictly */
		trans_new[iv] = dt*(channel_melt + channel_cavity + channel_creep) + trans;
		if (trans_new[iv] > Tmax)
		 trans_new[iv] = Tmax;
		if (trans_new[iv] < Tmin)
		 trans_new[iv] = Tmin;
	}

	/* Assign value */
	element->AddBasalInput(HydrologyTransmissivityEnum,trans_new,FINITEELEMENT);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(trans_new);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/

void HydrologyCuasAnalysis::UpdateEffectiveAquiferProperties(FemModel* femmodel){ /*{{{ */
	for(Object* & object : femmodel->elements->objects){
		Element *element = xDynamicCast<Element*>(object);
		UpdateEffectiveAquiferProperties(element);
	}
} /* }}} */
void HydrologyCuasAnalysis::UpdateEffectiveAquiferProperties(Element* element){ /* {{{ */

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries*/
	IssmDouble* Teff; // Effective transmissivity 
	IssmDouble* Seff; // Effective storage
	IssmDouble  Sprime;

	IssmDouble  psi;

	IssmDouble zb, head;
	IssmDouble conductivity;
	IssmDouble transmissivity;
	IssmDouble layer_thickness;
	IssmDouble ss, sy;

	int numnodes = basalelement->GetNumberOfNodes();

	Teff = xNew<IssmDouble>(numnodes);
	Seff = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	Input* base_input = basalelement->GetInput(BaseEnum); _assert_(base_input);
	Input* head_input = basalelement->GetInput(HydrologyHeadEnum); _assert_(head_input);
	Input* layer_thickness_input = basalelement->GetInput(HydrologySheetThicknessEnum); _assert_(layer_thickness_input);
	Input* trans_input= basalelement->GetInput(HydrologyTransmissivityEnum); _assert_(trans_input);
	Input* conductivity_input = basalelement->GetInput(HydrologyConductivityEnum); _assert_(conductivity_input);

	Input* ss_input = basalelement->GetInput(HydrologySpecificStorageEnum); _assert_(ss_input);
	Input* sy_input = basalelement->GetInput(HydrologySpecificYieldEnum); _assert_(sy_input);

	bool isconfined;
	IssmDouble  unconfinedSmooth;
	basalelement->FindParam(&isconfined,HydrologyIsConfinedEnum);
	basalelement->FindParam(&unconfinedSmooth,HydrologyCuasUnconfinedSmoothEnum);

	/*Start looping */
	Gauss *gauss=basalelement->NewGauss();
	for(int iv=0;iv<numnodes;iv++){

		gauss->GaussNode(basalelement->GetElementType(),iv);
		base_input->GetInputValue(&zb,gauss);
		head_input->GetInputValue(&head,gauss);
		trans_input->GetInputValue(&transmissivity,gauss);
		conductivity_input->GetInputValue(&conductivity,gauss);
		layer_thickness_input->GetInputValue(&layer_thickness,gauss);

		ss_input->GetInputValue(&ss,gauss);
		sy_input->GetInputValue(&sy,gauss);
		
		psi = head - zb;
		/*Neglect negative values*/
		if (psi < 0.0) psi = 0.01;

		/*First, guess effective storage coefficient.*/
		Seff[iv] = ss*layer_thickness;
		if (!isconfined){
			if (unconfinedSmooth > 0.0){
				if(layer_thickness <= psi){
				 Sprime = 0.0;
				}
				else if( (layer_thickness-unconfinedSmooth <= psi) & (psi<layer_thickness)){
					Sprime = sy/unconfinedSmooth*(layer_thickness-psi);
				}else{
					Sprime = sy;
				}
			} else {
				if (layer_thickness <= psi) Sprime = 0.0;
				else Sprime = sy;
			}
			Seff[iv] += Sprime;
		}

		/*Layer thickness from transmissivity and conductivity*/
		if (isconfined){
			Teff[iv]= transmissivity;
		}else{
			if (layer_thickness <= psi){
				Teff[iv] = transmissivity;
			}
			else{
				Teff[iv] = transmissivity * (psi/layer_thickness);
			}
		}

		if(xIsInf<IssmDouble>(Teff[iv])) _error_("Found inf value in Teff: " << Teff[iv] << "\n");
		if(xIsInf<IssmDouble>(Seff[iv])) _error_("Found inf value in Seff: " << Seff[iv] << "\n");
	}

	element->AddBasalInput(HydrologyTransmissivityEffectiveEnum,Teff,FINITEELEMENT); //element->GetElementType());
	element->AddBasalInput(HydrologyStorageEnum,Seff,FINITEELEMENT); //element->GetElementType());

	xDelete<IssmDouble>(Teff);
	xDelete<IssmDouble>(Seff);
	if(element->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
} /* }}} */

void HydrologyCuasAnalysis::UpdateEffectivePressure(FemModel* femmodel){ /*{{{*/
	for(Object* & object : femmodel->elements->objects){
		Element *element = xDynamicCast<Element*>(object);
		UpdateEffectivePressure(element);
	}
} /*}}}*/
void HydrologyCuasAnalysis::UpdateEffectivePressure(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;

	/*Intermediaries*/
	IssmDouble bed,thickness,head,layer_thickness;

	/* Fetch number of nodes and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* N = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	IssmDouble  g			 = element->FindParam(ConstantsGEnum);
	IssmDouble  rho_ice	 = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water  = element->FindParam(MaterialsRhoFreshwaterEnum);
	Input* head_input		= element->GetInput(HydrologyHeadEnum); _assert_(head_input);
	Input* thickness_input = element->GetInput(ThicknessEnum);	  _assert_(thickness_input);
	Input* base_input		= element->GetInput(BaseEnum);			 _assert_(base_input);
	Input* layer_thickness_input=element->GetInput(HydrologySheetThicknessEnum); _assert_(layer_thickness_input);

	Gauss* gauss=element->NewGauss();
	for (int iv=0;iv<numnodes;iv++){
		gauss->GaussNode(element->GetElementType(),iv);

		base_input->GetInputValue(&bed,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		head_input->GetInputValue(&head,gauss);
		layer_thickness_input->GetInputValue(&layer_thickness,gauss);

		N[iv] = rho_ice*g*thickness - rho_water*g*(head-bed-layer_thickness);
	}

	/*Set to 0 if inactive element*/
	if(element->IsAllFloating() || !element->IsIceInElement()){
		for(int iv=0;iv<numnodes;iv++) N[iv] = 0.;
		element->AddInput(EffectivePressureEnum,N,FINITEELEMENT);
		xDelete<IssmDouble>(N);
		return;
	}

	/*Add new gap as an input*/
	element->AddBasalInput(EffectivePressureEnum,N,FINITEELEMENT);

	/*Clean up and return*/
	xDelete<IssmDouble>(N);
	delete gauss;
}/*}}}*/

void HydrologyCuasAnalysis::ComputeWaterflux(FemModel* femmodel){/*{{{*/
	for(Object* & object : femmodel->elements->objects){
		Element *element = xDynamicCast<Element*>(object);
		ComputeWaterflux(element);
	}
}/*}}}*/
void HydrologyCuasAnalysis::ComputeWaterflux(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries*/
	IssmDouble* xyz_list=NULL;
	IssmDouble  Jdet;
	IssmDouble  dh[2];
	IssmDouble  dhdx, dhdy;
	IssmDouble  transmissivity;

	/* Fetch number of nodes and allocate output*/
	int numnodes = basalelement->GetNumberOfNodes();
	IssmDouble *flux, *fluxX, *fluxY;
	flux  = xNew<IssmDouble>(numnodes);
	fluxX = xNew<IssmDouble>(numnodes);
	fluxY = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	Input* head_input  = basalelement->GetInput(HydrologyHeadEnum); _assert_(head_input);
	Input* trans_input = basalelement->GetInput(HydrologyTransmissivityEffectiveEnum); _assert_(trans_input);

	Gauss* gauss=basalelement->NewGauss();
	for (int iv=0;iv<numnodes;iv++){
		gauss->GaussNode(basalelement->GetElementType(),iv);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		trans_input->GetInputValue(&transmissivity,gauss);
		head_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
		dhdx=dh[0];
		dhdy=dh[1];
		
		/*Now, compute water flux in x, y direction*/
		fluxX[iv] = -transmissivity*dhdx;
		fluxY[iv] = -transmissivity*dhdy;
		flux[iv]  = transmissivity*sqrt(dhdx*dhdx + dhdy*dhdy);
	}

	/*Add new gap as an input*/
	element->AddBasalInput(HydrologyBasalFluxEnum,flux,FINITEELEMENT);
	element->AddBasalInput(HydrologyWaterVxEnum,fluxX,FINITEELEMENT);
	element->AddBasalInput(HydrologyWaterVyEnum,fluxY,FINITEELEMENT);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(fluxX);
	xDelete<IssmDouble>(fluxY);
	xDelete<IssmDouble>(flux);
	delete gauss;
	if(element->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};

}/*}}}*/
