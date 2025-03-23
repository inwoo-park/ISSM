#include "./BasalforcingsLaddieSaltAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/TransientInput.h"

#define FINITEELEMENT P1Enum

/*Model processing*/
void BasalforcingsLaddieSaltAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int stabilization;
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");

	/*Do not add constraints in DG,  they are weakly imposed*/
	//if(stabilization!=3){
	//	IoModelToConstraintsx(constraints,iomodel,"md.basalforcings.spcthickness",MasstransportAnalysisEnum,FINITEELEMENT);
	//}
	
}/*}}}*/
void BasalforcingsLaddieSaltAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int penpair_ids[2];
	int count=0;
	int stabilization;
	int numvertex_pairing;

	/*Fetch parameters: */
	iomodel->FindConstant(&stabilization,"md.masstransport.stabilization");

	/*Loads only in DG*/
	if(stabilization==3){

		/*Get faces and elements*/
		CreateFaces(iomodel);
		iomodel->FetchData(1,"md.geometry.thickness");

		/*First load data:*/
		for(int i=0;i<iomodel->numberoffaces;i++){

			/*Get left and right elements*/
			int element=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]

			/*Now, if this element is not in the partition, pass: */
			if(!iomodel->my_elements[element]) continue;
			loads->AddObject(new Numericalflux(i+1,i,i,iomodel));
		}

		/*Free data: */
		iomodel->DeleteData(1,"md.geometry.thickness");
	}

	/*Create Penpair for vertex_pairing: */
	IssmDouble *vertex_pairing=NULL;
	IssmDouble *nodeonbase=NULL;
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(&nodeonbase,NULL,NULL,"md.mesh.vertexonbase");

	for(int i=0;i<numvertex_pairing;i++){

		if(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+0])-1]){

			/*In debugging mode, check that the second node is in the same cpu*/
			_assert_(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+1])-1]);

			/*Skip if one of the two is not on the bed*/
			if(iomodel->domaintype!=Domain2DhorizontalEnum){
				if(!(reCast<bool>(nodeonbase[reCast<int>(vertex_pairing[2*i+0])-1])) || !(reCast<bool>(nodeonbase[reCast<int>(vertex_pairing[2*i+1])-1]))) continue;
			}

			/*Get node ids*/
			penpair_ids[0]=reCast<int>(vertex_pairing[2*i+0]);
			penpair_ids[1]=reCast<int>(vertex_pairing[2*i+1]);

			/*Create Load*/
			loads->AddObject(new Penpair(count+1,&penpair_ids[0]));
			count++;
		}
	}

	/*Free resources: */
	iomodel->DeleteData(vertex_pairing,"md.masstransport.vertex_pairing");
	iomodel->DeleteData(nodeonbase,"md.mesh.vertexonbase");
}/*}}}*/
void BasalforcingsLaddieSaltAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  stabilization;
	iomodel->FindConstant(&stabilization,"md.masstransport.stabilization");

	/*Check in 3d*/
	if(stabilization==3 && iomodel->domaintype==Domain3DEnum) _error_("DG 3d not implemented yet");

	/*Create Nodes either DG or CG depending on stabilization*/
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	if(stabilization!=3){
		::CreateNodes(nodes,iomodel,BasalforcingsLaddieSaltAnalysisEnum,FINITEELEMENT,isamr);
	}
	else{
		::CreateNodes(nodes,iomodel,BasalforcingsLaddieSaltAnalysisEnum,P1DGEnum,isamr);
	}
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  BasalforcingsLaddieSaltAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void BasalforcingsLaddieSaltAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*
	NOTE: Nothing to do. 
	 */
}/*}}}*/
void BasalforcingsLaddieSaltAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	/*
	NOTE: Nothing to do. 
	 */
}/*}}}*/

/*Finite Element Analysis*/
void           BasalforcingsLaddieSaltAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BasalforcingsLaddieSaltAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BasalforcingsLaddieSaltAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BasalforcingsLaddieSaltAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BasalforcingsLaddieSaltAnalysis::CreateKMatrix(Element* element){/*{{{*/

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
ElementVector* BasalforcingsLaddieSaltAnalysis::CreatePVector(Element* element){/*{{{*/

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
void           BasalforcingsLaddieSaltAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,BasalforcingsLaddieTEnum);
}/*}}}*/
void           BasalforcingsLaddieSaltAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           BasalforcingsLaddieSaltAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Only update if on base*/
	if(!element->IsOnBase()) return;

	/*deal with logic of accumulating thickness if we are coupled to a 
	 * sea level core:*/
	int frequency,count,isgrd;
	element->FindParam(&isgrd,SolidearthSettingsGRDEnum);
	if(isgrd){
		element->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);
		element->FindParam(&count,SealevelchangeRunCountEnum);
	}

	/*Fetch dof list and allocate solution vector*/
	int *doflist = NULL;
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);

	int numnodes = element->GetNumberOfNodes();
	Input *thickness_input = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thickness_input);
	IssmDouble* Snew = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		Snew[i]=solution[doflist[i]];
		/*Check solution*/
		if(xIsNan<IssmDouble>(Snew[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(Snew[i])) _error_("Inf found in solution vector");
	}
	element->AddBasalInput(BasalforcingsLaddieTEnum,Snew,element->GetElementType());

	xDelete<int>(doflist);
	xDelete<IssmDouble>(Snew);
}/*}}}*/
void           BasalforcingsLaddieSaltAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/

/*Salt tranport analysis*/
ElementMatrix* BasalforcingsLaddieSaltAnalysis::CreateKMatrixCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement()) return NULL;

	/*Intermediaries */
	int        stabilization;
	int        domaintype,dim;
	IssmDouble Jdet,D_scalar,dt,h,factor;
	IssmDouble thickness; /*plume thickness*/
	IssmDouble vel,vx,vy,dvxdx,dvydy;
	IssmDouble xi,tau;
	IssmDouble dvx[2],dvy[2];
	IssmDouble D[4];
	IssmDouble* xyz_list = NULL;

	IssmDouble Kh; /*Horizontal diffusion coefficient*/

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		//case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		//case Domain3DEnum:           dim = 2; break;
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
	element->FindParam(&stabilization,MasstransportStabilizationEnum);

	element->FindParam(&Kh,BasalforcingsLaddieHorizontalDiffusivityEnum);

	/*Plume thickness D and depth averaged horizontal velocity*/
	Input* thickness_input = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thickness_input);
	Input* S_input  = element->GetInput(BasalforcingsLaddieSEnum); _assert_(S_input);
	Input* vx_input = element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
	Input* vy_input = element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);

	h = element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		/*Prepare inputs: */
		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);

		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

		thickness_input->GetInputValue(&thickness,gauss);
		dvxdx=dvx[0];
		dvydy=dvy[1];

		/*Transient term*/
		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*thickness*basis[i]*basis[j];

		/*Diffusion term: */
		factor = D_scalar*Kh*thickness;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += factor*(
							dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+j]
							);
			}
		}


		/*Advection term: */
		D_scalar=dt*gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/*\phi_i \phi_j \nabla\cdot v*/
				Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy);
				/*\phi_i v\cdot\nabla\phi_j*/
				Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*thickness*dbasis[0*numnodes+j] + vy*thickness*dbasis[1*numnodes+j]);
			}
		}
	}

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
			if(dim==1){
				vel=fabs(vx)+1.e-8;
			}
			else{
				vy_input->GetInputAverage(&vy);
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
			}
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

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* BasalforcingsLaddieSaltAnalysis::CreatePVectorCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int			stabilization,dim,domaintype;
	int         melt_style,point1;
	int         isentrainment;
	bool        mainlyfloating;
	IssmDouble  fraction1,fraction2;
	IssmDouble  Jdet,dt,intrusiondist;
	IssmDouble  thickness; /*plume thickness (or depth) [m]*/
	IssmDouble  S; /*plume salinity*/
	IssmDouble  Sa; /*ambient salinity*/
	IssmDouble  gammaS; /*salt exchange coefficient*/
	IssmDouble  vx,vy,vel,dvxdx,dvydy,xi,h,tau;
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
		//case Domain3DEnum:           dim = 2; break;
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
	//element->FindParam(&melt_style,GroundinglineMeltInterpolationEnum);
	element->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationEnum);

	Input* thickness_input= element->GetInput(BasalforcingsLaddieThicknessEnum);   _assert_(thickness_input);
	Input* S_input        = element->GetInput(BasalforcingsLaddieSEnum); _assert_(S_input);
	Input* entr_input     = element->GetInput(BasalforcingsLaddieEntrainmentRateEnum); _assert_(melt_input);
	/*Ambient temperature*/
	Input* Sa_input       = element->GetInput(BasalforcingsLaddieAmbientSalinityEnum); _assert_(Sa_input);
	h=element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/*Get input values*/
		thickness_input->GetInputValue(&thickness,gauss); /*plumem thickness*/
		S_input->GetInputValue(&S,gauss); /*plume salinity*/
		Sa_input->GetInputValue(&Sa,gauss); /*ambient salinity*/
		entr_input->GetInputValue(&entr_rate,gauss);

		IssmDouble factor = Jdet*gauss->weight*(S+dt*(entr_rate*Sa));
		for(int i=0;i<numnodes;i++) pe->values[i]+=factor*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return pe;
}/*}}}*/
