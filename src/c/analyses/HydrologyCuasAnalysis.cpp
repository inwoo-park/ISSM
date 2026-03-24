#include "./HydrologyCuasAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyCuasAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologycuasEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spchead",HydrologyCuasAnalysisEnum,P1Enum);
	
}/*}}}*/
void HydrologyCuasAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	return;

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

	/*Add input to elements*/
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	
    iomodel->FetchDataToInput(inputs,elements,"md.hydrology.conductivity",HydrologySheetConductivityEnum);
    iomodel->FetchDataToInput(inputs,elements,"md.hydrology.layer_thickness",HydrologySheetThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.Tinit",HydrologyTransmissivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.ss",HydrologySpecificStorageEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sy",HydrologySpecificYieldEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_heights",HydrologyBumpHeightEnum);
    iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_spacing",HydrologyBumpSpacingEnum);

    iomodel->FetchDataToInput(inputs,elements,"md.hydrology.moulin_input",HydrologyMoulinInputEnum);
    iomodel->FetchDataToInput(inputs,elements,"md.hydrology.englacial_input",HydrologyEnglacialInputEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vel",VelEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);

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
    Input* Teff_input    = basalelement->GetInput(HydrologyTransmissivityEffectiveEnum); _assert_(Teff_input);

	/* Start  looping on the number of gaussian points: */
    Gauss* gauss=basalelement->NewGauss(2);
    while(gauss->next()){

        basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
        basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
        basalelement->NodalFunctions(basis,gauss);

        storage_input->GetInputValue(&storage,gauss);
        Teff_input->GetInputValue(&Teff,gauss);

        /* Diffusion term */
        IssmDouble factor1 = dt*Teff*gauss->weight*Jdet;
        
        /* Transient mass term */
        IssmDouble factor2 = storage*gauss->weight*Jdet;

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
    IssmDouble* meltrate;
    IssmDouble* englacial_rate;

    /*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

    /*Initialize Element vector and other vectors*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

    /*Retrieve all inputs and parameters*/
    Input* head_input      = basalelement->GetInput(HydrologyHeadEnum); _assert_(head_input);
    Input* meltrate_input = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(meltrate_input);
    Input* englacial_input = basalelement->GetInput(HydrologyEnglacialInputEnum); _assert_(englacial_input);

    /*Get Params*/
    IssmDouble dt;
    basalelement->FindParam(&dt,TimesteppingTimeStepEnum);

    /* Start  looping on the number of gaussian points: */
    Gauss* gauss=basalelement->NewGauss(2);
    while(gauss->next()){

        basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

        meltrate_input->GetInputValue(&meltrate,gauss);
        englacial_input->GetInputValue(&englacial_rate,gauss);

        IssmDouble factor = (meltrate + englacial_rate)*Jdet*gauss->weight;

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
void HydrologyCuasAnalysis::UpdateTransmissivity(FemModel* femmodel){/*{{{*/
    for(Object* & object : femmodel->elements->objects){
        Element *element = xDynamicCast<Element*>(object);
        UpdateTransmissivity(element);
    }
}/*}}}*/
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
    IssmDouble  Jdet;

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
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  rho_ice   = basalelement->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = basalelement->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble  latent    = basalelement->FindParam(MaterialsLatentheatEnum);

	const int  numvertices= basalelement->GetNumberOfVertices();
	IssmDouble* trans_new = xNew<IssmDouble>(numvertices);

	IssmDouble vel;
	Input *vel_input  = basalelement->GetInput(VelEnum); _assert_(vel_input);
	Input *head_input = basalelement->GetInput(HydrologyHeadEnum);      _assert_(head_input);
	Input *B_input    = basalelement->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	Input *n_input    = basalelement->GetInput(MaterialsRheologyNEnum); _assert_(n_input);

	Input *conductivity_input=basalelement->GetInput(HydrologyConductivityEnum); _assert_(conductivity_input);
	Input *trans_input = basalelement->GetInput(HydrologyTransmissivityEnum); _assert_(trans_input);
	Input *Pe_input    = basalelement->GetInput(EffectivePressureEnum); _assert_(Pe_input);

	Input *bump_height_input = basalelement->GetInput(HydrologyBumpHeightEnum); _assert_(bump_height_input);
	Input *bump_space_input  = basalelement->GetInput(HydrologyBumpSpacingEnum); _assert_(bump_space_input);

	/* Update transmissivity */
	Gauss *gauss=basalelement->NewGauss(1);
    for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		trans_input->GetInputValue(&trans,gauss);
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
void HydrologyCuasAnalysis::UpdateStorage(FemModel* femmodel){/*{{{*/
    for(Object* & object : femmodel->elements->objects){
        Element *element = xDynamicCast<Element*>(object);
        UpdateStorage(element);
    }
}/*}}}*/
void HydrologyCuasAnalysis::UpdateStorage(Element* element){ /*{{{ */
	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	base_input=basalelement->GetInput(BaseEnum); _assert_(base_input);

	element->AddBasalInput(HydrologyStorageEnum,storage,P1DGEnum);
} /* }}} */
void HydrologyCuasAnalysis::UpdateEffectiveAquiferProperties(Element *element){ /* {{{ */

    /*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

    /*Intermediaries*/
    IssmDouble* Teff; // Effective transmissivity 
    IssmDouble* Seff; // Effective storage

    IssmDouble  b; // Layer thickness [m]
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
    Input* conductivity_input = basalelement->GetInput(HydrologySheetConductivityEnum); _assert_(conductivity_input);

    Input* ss_input = basalelement->GetInput(HydrologySpecificStorageEnum); _assert_(ss_input);
    Input* sy_input = basalelement->GetInput(HydrologySpecificYieldEnum); _assert_(sy_input);

    /*Start looping */
    Gauss *gauss=basalelement->NewGauss();
    for(int iv=0;iv<numnodes;iv++){

        gauss->GaussVertex(iv);
        base_input->GetInputValue(&zb,gauss);
        head_input->GetInputValue(&head,gauss);
        trans_input->GetInputValue(&transmissivity,gauss);
        conductivity_input->GetInputValue(&conductivity,gauss);

        layer_thickness_input->GetInputValue(&layer_thickness,gauss);

        ss_input->GetInputValue(&ss,gauss);
        sy_input->GetInputValue(&sy,gauss);
        
        psi = head - zb;
        /*Neglect negative values*/
        if (psi < 0.0) psi = 0.0;

        /*First, guess effective storage coefficient.*/
        Seff[iv] = ss + sy*layer_thickness;

        /*Layer thickness from transmissivity and conductivity*/
        // b = transmissivity/conductivity;
        if (psi < layer_thickness){
            Teff[iv] = transmissivity * (psi/layer_thickness);
        }
        else{
            Teff[iv] = transmissivity*psi/b;
        }
    }

    element->AddBasalInput(HydrologyTransmissivityEffectiveEnum,Teff,element->GetElementType());

    xDelete<IssmDouble>(Teff);
    xDelete<IssmDouble>(Seff);
    if(element->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
} /* }}} */
