#include "./BasalforcingsLaddieMassAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/TransientInput.h"

#define FINITEELEMENT P1Enum

/*Model processing*/
void BasalforcingsLaddieMassAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*No constraints for now*/
}/*}}}*/
void BasalforcingsLaddieMassAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads*/
}/*}}}*/
void BasalforcingsLaddieMassAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  stabilization;
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");

	/*Check in 3d*/
	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,BasalforcingsLaddieMassAnalysisEnum,FINITEELEMENT,isamr);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  BasalforcingsLaddieMassAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void BasalforcingsLaddieMassAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    isgroundingline;
	int    stabilization;
	int    finiteelement;
	int    grdmodel;

	/*Fetch data needed: */
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");

	/*Finite element type*/
	finiteelement = FINITEELEMENT;
	if(stabilization==3){
		finiteelement = P1DGEnum;
	}

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	/*Create Inptus: */
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.D",BasalforcingsLaddieThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.vx",BasalforcingsLaddieVyEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.vy",BasalforcingsLaddieVxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.temperature",BasalforcingsLaddieTEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.salinity",BasalforcingsLaddieSEnum);

	if(isgroundingline) 	iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);

	/*Initialize LaddieThicknessResidualEnum input*/
	InputUpdateFromConstantx(inputs,elements,0.,BasalforcingsLaddieThicknessResidualEnum);

	/*Get what we need for ocean-induced basal melting*/
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}

	/*Deal with forcing fields: temperature and salinity*/
	IssmDouble* array2d_T=NULL;
	IssmDouble* array2d_S=NULL;
	IssmDouble* temp=NULL;
	int M,N;
	int num_depths;

	iomodel->FetchData(&temp,&M,&num_depths,"md.basalforcings.forcing_depth");
	xDelete<IssmDouble>(temp);
	_assert_(M==1); _assert_(num_depths>=1);

	for (int kk=0;kk<num_depths;kk++){

		/*Fetch forcing temperature and salinity*/
		iomodel->FetchData(&array2d_T, &M, &N, kk, "md.basalforcings.forcing_temperature");
		if(!array2d_T) _error_("md.basalforcings.forcing_temperature not found in binary file");
		iomodel->FetchData(&array2d_S, &M, &N, kk, "md.basalforcings.forcing_salinity");
		if(!array2d_S) _error_("md.basalforcings.forcing_salinity not found in binary file");

		for(Object* & object : elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			if(iomodel->domaintype!=Domain2DhorizontalEnum && !element->IsOnBase()) continue;
			element->DatasetInputAdd(BasalforcingsLaddieForcingTemperatureEnum,array2d_T,inputs,iomodel,M,N,1,BasalforcingsLaddieForcingTemperatureEnum,kk);
			element->DatasetInputAdd(BasalforcingsLaddieForcingSalinityEnum,array2d_T,inputs,iomodel,M,N,1,BasalforcingsLaddieForcingSalinityEnum,kk);
		}
	}
	xDelete<IssmDouble>(array2d_T);
	xDelete<IssmDouble>(array2d_S);
}/*}}}*/
void BasalforcingsLaddieMassAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	//int     numoutputs;
	//char**  requestedoutputs = NULL;
	IssmDouble *transparam;
	int M, N;

	/*
	!!NOTE: Parameters for LADDIE are updated through "CreateParameters.cpp"
	Skip this part?!
	Recheck this part in further.
	 */

	/*Basalforcing plume model specials*/
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Kh",BasalforcingsLaddieHorizontalDiffusivityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Ah",BasalforcingsLaddieHorizontalViscosityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Dmin",BasalforcingsLaddieThicknessMinEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Utide",BasalforcingsLaddieVelTideEnum));

	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Cd",BasalforcingsLaddieCdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Cd_top",BasalforcingsLaddieCdTopEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.f_cori",BasalforcingsLaddieCoriolisFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Kparam",BasalforcingsLaddieKparamEnum));

	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.ismelt",BasalforcingsLaddieIsMeltEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.isentrainment",BasalforcingsLaddieIsEntrainmentEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.stabilization",BasalforcingsLaddieStabilizationEnum));

	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.subtimestep",BasalforcingsLaddieSubTimestepEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.subtimestep",BasalforcingsLaddieSubTimestepDummyEnum));

	/*Deal with forcing depth: */
	iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.forcing_depth");
	_assert_(M==1); _assert_(N>=1);
	parameters->AddObject(new DoubleVecParam(BasalforcingsLaddieForcingDepthEnum,transparam,N));
	xDelete<IssmDouble>(transparam);

	//iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.basalforcings.requested_outputs");
	//parameters->AddObject(new IntParam(MasstransportNumRequestedOutputsEnum,numoutputs));
	//if(numoutputs)parameters->AddObject(new StringArrayParam(MasstransportRequestedOutputsEnum,requestedoutputs,numoutputs));
	//iomodel->DeleteData(&requestedoutputs,numoutputs,"md.basalforcings.requested_outputs");
}/*}}}*/

/*Finite Element Analysis*/
void           BasalforcingsLaddieMassAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BasalforcingsLaddieMassAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BasalforcingsLaddieMassAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BasalforcingsLaddieMassAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BasalforcingsLaddieMassAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(!element->IsOceanInElement()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementMatrix* Ke = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			Ke = CreateKMatrixCG(basalelement);
			break;
		case P0DGEnum:
		case P1DGEnum:
			//Ke = CreateKMatrixDG(basalelement);
			_error_("Not implemented yet.");
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};

	return Ke;
}/*}}}*/
ElementVector* BasalforcingsLaddieMassAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementVector* pe = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			pe = CreatePVectorCG(basalelement);
			break;
		case P0DGEnum:
		case P1DGEnum:
			//pe = CreatePVectorDG(basalelement);
			_error_("Not implemented yet");
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void           BasalforcingsLaddieMassAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,BasalforcingsLaddieThicknessEnum);
}/*}}}*/
void           BasalforcingsLaddieMassAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           BasalforcingsLaddieMassAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Only update if on base*/
	if(!element->IsOnBase() || !element->IsIceInElement() || !element->IsOceanInElement()) return;

	/*Fetch dof list and allocate solution vector*/
	int *doflist = NULL;
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);

	int numnodes = element->GetNumberOfNodes();
	IssmDouble* Dnew      = xNew<IssmDouble>(numnodes);
	IssmDouble* Dresidual = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	IssmDouble Dmin = element->FindParam(BasalforcingsLaddieThicknessMinEnum);
	for(int i=0;i<numnodes;i++){
		Dnew[i]=solution[doflist[i]];
		Dresidual[i]=0.;
		/*Check solution*/
		if(xIsNan<IssmDouble>(Dnew[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(Dnew[i])) _error_("Inf found in solution vector");
		if(Dnew[i]<Dmin){
			Dresidual[i]=Dmin-Dnew[i];
			Dnew[i]=Dmin;
		}
	}
	element->AddBasalInput(BasalforcingsLaddieThicknessEnum,Dnew,element->GetElementType());
	element->AddBasalInput(BasalforcingsLaddieThicknessResidualEnum,Dresidual,element->GetElementType());

	xDelete<int>(doflist);
	xDelete<IssmDouble>(Dnew);
 	xDelete<IssmDouble>(Dresidual);
}/*}}}*/
void           BasalforcingsLaddieMassAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Deal with ocean constraint*/
	SetActiveNodesLSMx(femmodel);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);
		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

		IssmDouble Dmin = element->FindParam(BasalforcingsLaddieThicknessMinEnum);
		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&ls_active[0],IceMaskNodeActivationEnum);

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]<0. && ls_active[in]==1.){
				node->Activate();
			}
			else{
				/*Apply plume thickness as zero value along grounding line*/
				node->Deactivate();
				node->ApplyConstraint(0,Dmin);
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(ls_active);
	}

	return;
}/*}}}*/

/*Mass balance special*/
ElementMatrix* BasalforcingsLaddieMassAnalysis::CreateKMatrixCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int        stabilization;
	int        domaintype,dim;
	IssmDouble Jdet,D_scalar,dt,h,factor;
	IssmDouble vel,vx,vy,dvxdx,dvydy;
	IssmDouble melt, entrain; /*unit m s-1*/
	IssmDouble xi,tau;
	IssmDouble dvx[2],dvy[2];
	IssmDouble D[4];
	IssmDouble* xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		//case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*		dbasis = xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationEnum);
	/*Plume thickness D and depth averaged horizontal velocity*/
	Input* vx_input  = element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
	Input* vy_input  = element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);

	h = element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		/*Transient term*/
		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j];

		/*Advection terms*/
		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);

		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

		D_scalar=dt*gauss->weight*Jdet;
		dvxdx=dvx[0];
		dvydy=dvy[1];
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/*\phi_i \phi_j \nabla\cdot v*/
				Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy);
				/*\phi_i v\cdot\nabla\phi_j*/
				Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
			}
		}

		/*Prepare parameters for advection stabilization scheme*/
		for(int i=0;i<4;i++) D[i]=0.;
		switch(stabilization){
			case 0:
				/*Nothing to be done*/
				break;
			case 1:
				/*SSA*/
				vx_input->GetInputAverage(&vx);
				if(dim==2) vy_input->GetInputAverage(&vy);
				D[0*dim+0]=h/2.0*fabs(vx);
				if(dim==2) D[1*dim+1]=h/2.0*fabs(vy);
				break;
			case 2:
				/*Streamline upwinding*/
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				tau=h/(2*vel);
				break;
			case 5:
				/*SUPG*/
				if(dim!=2) _error_("Stabilization "<<stabilization<<" not supported yet for dim != 2");
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				//xi=0.3130;
				xi=1;
				tau=xi*h/(2*vel);
				//tau=dt/6; // as implemented in Ua
				break;
			default:
				_error_("Stabilization "<<stabilization<<" not supported yet");
		}

		if(stabilization==1){/*{{{*/
			/*SSA*/
			if(dim==1) D[0]=D_scalar*D[0];
			else{
				D[0*dim+0]=D_scalar*D[0*dim+0];
				D[1*dim+0]=D_scalar*D[1*dim+0];
				D[0*dim+1]=D_scalar*D[0*dim+1];
				D[1*dim+1]=D_scalar*D[1*dim+1];
			}

			if(dim==2){
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += (
									dbasis[0*numnodes+i] *(D[0*dim+0]*dbasis[0*numnodes+j] + D[0*dim+1]*dbasis[1*numnodes+j]) +
									dbasis[1*numnodes+i] *(D[1*dim+0]*dbasis[0*numnodes+j] + D[1*dim+1]*dbasis[1*numnodes+j])
									);
					}
				}
			}
			else{
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += dbasis[0*numnodes+i]*D[0]*dbasis[0*numnodes+j];
					}
				}
			}
		}/*}}}*/
		if(stabilization==2){/*{{{*/
			/*Streamline upwind*/
			_assert_(dim==2);
			factor = dt*gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i])*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j]);
				}
			}
		}/*}}}*/
		if(stabilization==5){/*{{{*/
			/*Mass matrix - part 2*/
			factor = gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*basis[j]*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
				}
			}
			/*Mass matrix - part 3*/
			factor = gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*basis[j]*(basis[i]*dvxdx+basis[i]*dvydy);
				}
			}

			/*Advection matrix - part 2, A*/
			factor = dt*gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j])*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
				}
			}
			/*Advection matrix - part 3, A*/
			factor = dt*gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j])*(basis[i]*dvxdx+basis[i]*dvydy);
				}
			}

			/*Advection matrix - part 2, B*/
			factor = dt*gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(basis[j]*dvxdx+basis[j]*dvydy)*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
				}
			}
			/*Advection matrix - part 3, B*/
			factor = dt*gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(basis[j]*dvxdx+basis[j]*dvydy)*(basis[i]*dvxdx+basis[i]*dvydy);
				}
			}
		}/*}}}*/
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* BasalforcingsLaddieMassAnalysis::CreatePVectorCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int			stabilization,dim,domaintype;
	int         melt_style,point1;
	int         isentrainment;
	bool        mainlyfloating;
	IssmDouble  factor;
	IssmDouble  fraction1,fraction2;
	IssmDouble  Jdet,dt,intrusiondist;
	IssmDouble  ms,mb,gmb,fmb,fmb_pert,gldistance;
	IssmDouble  vx,vy,vel,dvxdx,dvydy,xi,h,tau;
	IssmDouble  thickness; /*Plume thickness (or depth) [m]*/
	IssmDouble  melt_rate, entr_rate; /*melting rate and entrainment rate [unit: m s-1]*/
	IssmDouble  g; /* gravitational acceleration [m s-1]*/
	IssmDouble  dvx[2],dvy[2];
	IssmDouble  gllevelset,phi=1.;
	IssmDouble* xyz_list = NULL;
	Gauss*      gauss     = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&isentrainment,BasalforcingsLaddieIsEntrainmentEnum);
	switch(domaintype){
		//case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis= xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationEnum);
	element->FindParam(&g,ConstantsGEnum);

	Input* thickness_input  = element->GetInput(BasalforcingsLaddieThicknessEnum);   _assert_(thickness_input);
	Input* vx_input      = element->GetInput(BasalforcingsLaddieVxEnum);		  _assert_(vx_input);
	Input* vy_input      = element->GetInput(BasalforcingsLaddieVyEnum);	     _assert_(vy_input);
	Input* melt_input    = element->GetInput(BasalforcingsFloatingiceMeltingRateEnum); _assert_(melt_input);
	Input* entr_input    = element->GetInput(BasalforcingsLaddieEntrainmentRateEnum); _assert_(entr_input);
	h=element->CharacteristicLength();

	/*Recover portion of element that is grounded*/
	phi=element->GetGroundedPortion(xyz_list);
	if(melt_style==SubelementMelt2Enum){
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,MaskOceanLevelsetEnum,0);
	    gauss = element->NewGauss(point1,fraction1,fraction2,3);
	}
	else if(melt_style==IntrusionMeltEnum){
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,DistanceToGroundinglineEnum,intrusiondist);
       	gauss = element->NewGauss(point1,fraction1,fraction2,3);
	}
	else{
		gauss = element->NewGauss(3);
	}

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/*Get input values*/
		thickness_input->GetInputValue(&thickness,gauss);
		melt_input->GetInputValue(&melt_rate,gauss);
		entr_input->GetInputValue(&entr_rate,gauss);

		factor = Jdet*gauss->weight*(thickness+dt*(melt_rate+entr_rate));
		for(int i=0;i<numnodes;i++) pe->values[i]+=factor*basis[i];

		if(stabilization==5){ //SUPG
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			vx_input->GetInputAverage(&vx);
			vy_input->GetInputAverage(&vy);
			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			vel=sqrt(vx*vx+vy*vy)+1.e-8;
			dvxdx=dvx[0];
			dvydy=dvy[1];
			//xi=0.3130;
			xi=1;
			tau=xi*h/(2*vel);
			//tau=dt/6; // as implemented in Ua

			/*Force vector - part 2*/
			factor = Jdet*gauss->weight*(thickness+dt*(melt_rate+entr_rate));
			for(int i=0;i<numnodes;i++){
				pe->values[i]+=factor*(tau*vx*dbasis[0*numnodes+i]+tau*vy*dbasis[1*numnodes+i]);
			}
			/*Force vector - part 3*/
			for(int i=0;i<numnodes;i++){
				pe->values[i]+=factor*(tau*basis[i]*dvxdx+tau*basis[i]*dvydy);
			}
		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return pe;
}/*}}}*/
