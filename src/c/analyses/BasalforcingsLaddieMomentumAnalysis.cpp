#include "./BasalforcingsLaddieMomentumAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/TransientInput.h"

#define FINITEELEMENT P1Enum

/*Model processing*/
void BasalforcingsLaddieMomentumAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*Nothing to do*/	
}/*}}}*/
void BasalforcingsLaddieMomentumAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int penpair_ids[2];
	int count=0;
	int stabilization;
	int numvertex_pairing;

	/*Fetch parameters: */
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilizationMomentum");

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
void BasalforcingsLaddieMomentumAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  stabilization;
	int  approximation;
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilizationMomentum");

	/*Check in 3d*/
	if(stabilization==3 && iomodel->domaintype==Domain3DEnum) _error_("DG 3d not implemented yet");

	/*Create Nodes either DG or CG depending on stabilization*/
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");

	/* Create 2 nodes for vx and vy*/
	int *approximations = xNew<int>(iomodel->numberofvertices);
	for(int i=0;i<iomodel->numberofvertices;i++){
	 approximations[i] = IoCodeToEnumVertexEquation(2);
	}

	if(stabilization!=3){
		::CreateNodes(nodes,iomodel,BasalforcingsLaddieMomentumAnalysisEnum,FINITEELEMENT,isamr,0,approximations);
	}
	else{
		::CreateNodes(nodes,iomodel,BasalforcingsLaddieMomentumAnalysisEnum,P1DGEnum,isamr,0,approximations);
	}
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");

	delete approximations;
}/*}}}*/
int  BasalforcingsLaddieMomentumAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	/*
	NOTE:This script follows the StressbalanceAnalysis::DofsPerNode

	Number of dofs for vx and vy = 2
	*/

	return 2;
}/*}}}*/
void BasalforcingsLaddieMomentumAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    stabilization;
	int    finiteelement;
	int   *finiteelement_list;

	/*Fetch data needed: */
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");
	finiteelement_list=xNewZeroInit<int>(iomodel->numberofelements);

	/*Finite element type*/
	//finiteelement = FINITEELEMENT;
	//if(stabilization==3){
	//	finiteelement = P1DGEnum;
	//}
	//for(int i=0;i<iomodel->numberofelements;i++){
	//	finiteelement_list[i]=finiteelement;
	//}

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Clear memory: */
	xDelete<int>(finiteelement_list);
}/*}}}*/
void BasalforcingsLaddieMomentumAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	/*
	NOTE: Nothing to do in this part. See "UpdateParameters" in "BasalforcingsLaddieMassAnalay.cpp"
	*/ 
}/*}}}*/

void           BasalforcingsLaddieMomentumAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BasalforcingsLaddieMomentumAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BasalforcingsLaddieMomentumAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BasalforcingsLaddieMomentumAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BasalforcingsLaddieMomentumAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement() || !element->IsOceanInElement()) return NULL;

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
ElementVector* BasalforcingsLaddieMomentumAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement() || !element->IsOceanInElement()) return NULL;

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
void           BasalforcingsLaddieMomentumAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	/*
	 Look up "GetSolutionFromInputsHoriz" in "StressbalanceAnalysis.cpp".
	 */

	IssmDouble vx, vy;
	int        domaintype, dim, approximation, dofpernode;
	int*       doflist = NULL;

	/*Get some parameters*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; dofpernode = 2; break;
		case Domain3DEnum:           dim = 3; dofpernode = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*dofpernode;

	/*Fetch dof list and allocate solution vector*/
	element->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numdof);

	/*Get inputs*/
	Input* vx_input=element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	Gauss* gauss=element->NewGauss();
	for(int i=0;i<numnodes;i++){
		gauss->GaussNode(element->FiniteElement(),i);

		/*Recover vx and vy*/
		vx_input->GetInputValue(&vx,gauss);
		values[i*dofpernode+0]=vx;

		vy_input->GetInputValue(&vy,gauss);
		values[i*dofpernode+1]=vy;
	}

	solution->SetValues(numdof,doflist,values,INS_VAL);

	/*Free resources:*/
	delete gauss;
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}/*}}}*/
void           BasalforcingsLaddieMomentumAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           BasalforcingsLaddieMomentumAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int         i,dim=2,domaintype;
	int*        doflist=NULL;
	IssmDouble* xyz_list=NULL;
	IssmDouble  vcut;
	Element*    basalelement=NULL;

	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&vcut, BasalforcingsLaddieVCutEnum);
	/*Get basal element*/
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!element->IsOnBase()){xDelete<IssmDouble>(xyz_list); return;}
			basalelement=element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();
	int numdof   = numnodes*dim;

	/*Fetch dof list and allocate solution vectors*/
	basalelement->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numdof);
	IssmDouble* vx     = xNew<IssmDouble>(numnodes);
	IssmDouble* vy     = xNew<IssmDouble>(numnodes);
	IssmDouble* vel    = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	//if(dim==2) basalelement->TransformSolutionCoord(&values[0],XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		vx[i]=values[i*dim+0];
		if(xIsNan<IssmDouble>(vx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i])) _error_("Inf found in solution vector");

		vy[i]=values[i*dim+1];
		if(xIsNan<IssmDouble>(vy[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vy[i])) _error_("Inf found in solution vector");
	}
	
	/*Update plume velocity:  */
	double epsilon=1e-6;
	for(i=0;i<numnodes;i++){
		vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
		/*Cut for stability*/
		if (vel[i] > vcut){
			vx[i] = vx[i]*vcut/vel[i];
			vy[i] = vy[i]*vcut/vel[i];
			vel[i] = vcut;
		}
	}

	/*Cut for stability*/
	//for(i=0;i<numnodes;i++){
	//	vx[i] = min(vx[i], vcut);
	//	vx[i] = max(vx[i], -vcut);

	//	vy[i] = min(vy[i], vcut);
	//	vy[i] = max(vy[i], -vcut);
	//}


	/*Add vx and vy as inputs to the tria element: */
	element->AddBasalInput(BasalforcingsLaddieVxEnum,vx,element->GetElementType());
	element->AddBasalInput(BasalforcingsLaddieVyEnum,vy,element->GetElementType());
	element->AddBasalInput(BasalforcingsLaddieVelEnum,vel,element->GetElementType());

	/*Free resources:*/
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           BasalforcingsLaddieMomentumAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	SetActiveNodesLSMx(femmodel);

	int         domaintype;
	Node       *node;
	Element    *element;
	int         numnodes;
	IssmDouble *mask;
	IssmDouble *ls_active;
	IssmDouble  spcvx=0.0, spcvy=0.0;

	IssmDouble *xyz_list=NULL;
	IssmDouble  xyz_front[2][3];
	IssmDouble *normal;

	/*Find domain type*/
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		element  = xDynamicCast<Element*>(object);
		numnodes = element->GetNumberOfNodes();

		/*Get xyz coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		mask      = xNew<IssmDouble>(numnodes);
		ls_active = xNew<IssmDouble>(numnodes);
		normal    = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&ls_active[0],IceMaskNodeActivationEnum);

		for(int in=0;in<numnodes;in++){
			node=element->GetNode(in);
			if(mask[in]<0. && ls_active[in]==1.){
				node->Activate();
			}
			else{
				/*Apply plume vector (vx, vy) as zero value along grounding line*/
				node->Deactivate();
				node->ApplyConstraint(0,spcvx); /*for vx*/
				node->ApplyConstraint(1,spcvy); /*for vy*/
			}
		}

		//if(false)/* NOTE: some contraints along the grounding line.. ?? {{{*/
		//	if(!element->IsIceInElement()){
		//		/*nothing to do*/
		//	}else if(!element->IsZeroLevelset(MaskOceanLevelsetEnum)){
		//		/*nothing to do*/
		//	}else{
		//		/*Search grounding line*/
		//		/*Be sure that values are not zero*/
		//		if(mask[0]==0.) mask[0]=mask[0]+epsilon;
		//		if(mask[1]==0.) mask[1]=mask[1]+epsilon;
		//		if(mask[2]==0.) mask[2]=mask[2]+epsilon;

		//		int pt1 = 0;
		//		int pt2 = 1;
		//		if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2

		//			/*Portion of the segments*/
		//			s1=gl[2]/(gl[2]-gl[1]);
		//			s2=gl[2]/(gl[2]-gl[0]);
		//			if(gl[2]<0.){
		//				pt1 = 1; pt2 = 0;
		//			}
		//			xyz_front[pt2][0]=xyz_list[2][0]+s1*(xyz_list[1][0]-xyz_list[2][0]);
		//			xyz_front[pt2][1]=xyz_list[2][1]+s1*(xyz_list[1][1]-xyz_list[2][1]);
		//			xyz_front[pt2][2]=xyz_list[2][2]+s1*(xyz_list[1][2]-xyz_list[2][2]);
		//			xyz_front[pt1][0]=xyz_list[2][0]+s2*(xyz_list[0][0]-xyz_list[2][0]);
		//			xyz_front[pt1][1]=xyz_list[2][1]+s2*(xyz_list[0][1]-xyz_list[2][1]);
		//			xyz_front[pt1][2]=xyz_list[2][2]+s2*(xyz_list[0][2]-xyz_list[2][2]);
		//		}
		//		else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2

		//			/*Portion of the segments*/
		//			s1=gl[0]/(gl[0]-gl[1]);
		//			s2=gl[0]/(gl[0]-gl[2]);
		//			if(gl[0]<0.){
		//				pt1 = 1; pt2 = 0;
		//			}

		//			xyz_front[pt1][0]=xyz_list[0][0]+s1*(xyz_list[1][0]-xyz_list[0][0]);
		//			xyz_front[pt1][1]=xyz_list[0][1]+s1*(xyz_list[1][1]-xyz_list[0][1]);
		//			xyz_front[pt1][2]=xyz_list[0][2]+s1*(xyz_list[1][2]-xyz_list[0][2]);
		//			xyz_front[pt2][0]=xyz_list[0][0]+s2*(xyz_list[2][0]-xyz_list[0][0]);
		//			xyz_front[pt2][1]=xyz_list[0][1]+s2*(xyz_list[2][1]-xyz_list[0][1]);
		//			xyz_front[pt2][2]=xyz_list[0][2]+s2*(xyz_list[2][2]-xyz_list[0][2]);
		//		}
		//		else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2

		//			/*Portion of the segments*/
		//			s1=gl[1]/(gl[1]-gl[0]);
		//			s2=gl[1]/(gl[1]-gl[2]);
		//			if(gl[1]<0.){
		//				pt1 = 1; pt2 = 0;
		//			}

		//			xyz_front[pt2][0]=xyz_list[1][0]+s1*(xyz_list[0][0]-xyz_list[1][0]);
		//			xyz_front[pt2][1]=xyz_list[1][1]+s1*(xyz_list[0][1]-xyz_list[1][1]);
		//			xyz_front[pt2][2]=xyz_list[1][2]+s1*(xyz_list[0][2]-xyz_list[1][2]);
		//			xyz_front[pt1][0]=xyz_list[1][0]+s2*(xyz_list[2][0]-xyz_list[1][0]);
		//			xyz_front[pt1][1]=xyz_list[1][1]+s2*(xyz_list[2][1]-xyz_list[1][1]);
		//			xyz_front[pt1][2]=xyz_list[1][2]+s2*(xyz_list[2][2]-xyz_list[1][2]);
		//		}
		//		else{
		//			_error_("case not possible");
		//		}

		//		/*Some checks in debugging mode*/
		//		_assert_(s1>=0 && s1<=1.);
		//		_assert_(s2>=0 && s2<=1.);

		//		/*Get normal vector*/
		//		element->NormalSection(&normal[0],&xyz_front[0][0]);
		//	}
		//}/*}}}*/

		/*Find two vertices pair*/	
		xDelete<IssmDouble>(normal);
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(ls_active);
		xDelete<IssmDouble>(xyz_list);
	}
}/*}}}*/

/*CG method*/
ElementMatrix* BasalforcingsLaddieMomentumAnalysis::CreateKMatrixCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         stabilization;
	int         domaintype,dim;
	IssmDouble  Jdet,D_scalar,h,factor;
	IssmDouble  dt; /*time-stepping*/
	IssmDouble  thickness;
	IssmDouble  vel,vx,vy,dvxdx,dvxdy,dvydx,dvydy;
	/*Hvel ~ thickness \times vel*/
	IssmDouble  vx_avg, vy_avg, thickness_avg;
	IssmDouble  Hvx, Hvy, Hvel;
	IssmDouble  melt, entrain; /*unit m s-1*/
	IssmDouble  Cd;
	IssmDouble  xi,tau;
	IssmDouble  dvx[2],dvy[2];
	IssmDouble  dthk[2];
	IssmDouble  dthkdx, dthkdy;
	IssmDouble  D[2][2]={0.0};
	IssmDouble *xyz_list = NULL;
	IssmDouble  g; /*gravitational acceleration [m s-2]*/

	IssmDouble  Ah; /*horizontal viscosity*/
	IssmDouble  coriolis_freq; /*Coriolis frequency [s-1]*/
	int         isconvection;

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
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationMomentumEnum);
	element->FindParam(&coriolis_freq,BasalforcingsLaddieCoriolisFrequencyEnum);
	element->FindParam(&Ah,BasalforcingsLaddieHorizontalViscosityEnum);
	element->FindParam(&Cd,BasalforcingsLaddieCdEnum);
	element->FindParam(&g,ConstantsGEnum);
	element->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);
	element->FindParam(&isconvection,BasalforcingsLaddieIsConvectionEnum);

	/*Plume thickness D and depth averaged horizontal velocity*/
	Input* thickness_input = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thickness_input);
	Input* vx_input        = element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
	Input* vy_input        = element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);

	h = element->CharacteristicLength();

	/* Start looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		
		/*Retrieve inputs: */
		thickness_input->GetInputValue(&thickness,gauss);
		thickness_input->GetInputAverage(&thickness_avg);
		thickness_input->GetInputDerivativeValue(&dthk[0],xyz_list,gauss);
		dthkdx = dthk[0];
		dthkdy = dthk[1];

		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputAverage(&vx_avg);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);

		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputAverage(&vy_avg);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

		dvxdx=dvx[0];
		dvxdy=dvx[1];
		dvydx=dvy[0];
		dvydy=dvy[1];

		vel = sqrt(pow(vx,2.0) + pow(vy,2.0));;

		/*Transient term: */
		factor=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[2*i*2*numnodes+2*j]       += factor*thickness*basis[i]*basis[j];
				Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*thickness*basis[i]*basis[j];
			}
		}

		/*Convection term: */
		D_scalar=gauss->weight*Jdet*dt;
		if(isconvection==0){
			/*NOTE: Nothing to be done*/
			//_printf0_("Momentum: Nothing to be done in convection term.\n");
		}
		else if(isconvection==1){
			/*
			NOTE: Get diagonal component
			| C_(2i-1, 2j-1) 0         |
			| 0              C_(2i,2j) |
			*/
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*term: dDu2/dx + dDuv/dy = 2Du dudx + u^2 dDdx + Du dvdy + Dv dudy + uv dDdy*/
					Ke->values[2*i*2*numnodes+2*j] += D_scalar*basis[i]*(
								+ 2*thickness*vx*dbasis[0*numnodes+j]
								+ vx*dthkdx*basis[j]
								+ thickness*dvydy*basis[j]
								+ thickness*vy*dbasis[1*numnodes+j]
								+ vy*dthkdy*basis[j]
								);

					/*term: dDuv/dx + dDv^2/dy = Du dvdx + Dv dudx + uv dDdx + 2Dv dvdy  + v^2 dDdy*/
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += D_scalar*basis[i]*(
								+thickness*vx*dbasis[0*numnodes+j]
								+thickness*dvxdx*basis[j]
								+vx*dthkdx*basis[j]
								+2*thickness*vy*dbasis[1*numnodes+j]
								+vy*dthkdy*basis[j]
								);

				}
			}

			//for(int i=0;i<numnodes;i++){
			//	for(int j=0;j<numnodes;j++){
			//		/*term: 2u dudx + u dvdy + v dudy*/
			//		Ke->values[2*i*2*numnodes+2*j] += D_scalar*thickness*basis[i]*basis[i]*(
			//					+ 2*vx*dbasis[0*numnodes+j] 
			//					+ vy*dbasis[1*numnodes+j]
			//					+ dvydy*basis[j]
			//					);

			//		/*term: u dvdx + v dudx + 2v dvdy -> u dvdx + 2v dvdy*/
			//		Ke->values[(2*i+1)*2*numnodes+2*j+1] += D_scalar*thickness*basis[i]*basis[i]*(
			//					+ dvxdx*basis[j]
			//					+ vx*dbasis[0*numnodes+j]
			//					+ 2*vy*dbasis[1*numnodes+j]
			//					);

			//	}
			//}
		}
		else if(isconvection==2){
			/*
			NOTE: 
			 */
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*term: 2u dudx + u dvdy + v dudy -> 2u dudx + v dudy*/
					Ke->values[2*i*2*numnodes+2*j] += D_scalar*(
								+ 2*vx*basis[i]*(thickness*dbasis[0*numnodes+j] + dthkdx*basis[j])
								+ vy*basis[i]*(thickness*dbasis[1*numnodes+j] + dthkdy*basis[j])
								);
					/*term: 2u dudx + u dvdy + v dudy -> u dvdy*/
					Ke->values[2*i*2*numnodes+2*j+1] += D_scalar*(
								basis[i]*vx*(thickness*dbasis[1*numnodes+j] + dthkdy*basis[j])
								);

					/*term: u dvdx + v dudx + 2v dvdy -> v dudx */
					Ke->values[(2*i+1)*2*numnodes+2*j] += D_scalar*(
								basis[i]*vy*(thickness*dbasis[0*numnodes+j] + dthkdx*basis[j])
								);
					/*term: u dvdx + v dudx + 2v dvdy -> u dvdx + 2v dvdy*/
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += D_scalar*(
								vx*basis[i]*(thickness*dbasis[0*numnodes+j] + dthkdx*basis[j])
								+ 2*basis[i]*vy*(thickness*dbasis[1*numnodes+j] + dthkdy*basis[j])
								);
				}
			}
		}
		else if(isconvection==3){
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*term: 2u dudx + u dvdy + v dudy -> 2u dudx + v dudy*/
					Ke->values[2*i*2*numnodes+2*j] += D_scalar*thickness*basis[i]*(
								+ 2*vx*basis[i]*dbasis[0*numnodes+j] 
								+ vy*basis[i]*dbasis[1*numnodes+j]
								);
					/*term: 2u dudx + u dvdy + v dudy -> u dvdy*/
					Ke->values[2*i*2*numnodes+2*j+1] += D_scalar*thickness*basis[i]*(
								vx*basis[i]*dbasis[1*numnodes+j] 
								);

					/*term: u dvdx + v dudx + 2v dvdy -> v dudx */
					Ke->values[(2*i+1)*2*numnodes+2*j] += D_scalar*thickness*basis[i]*(
								vy*basis[i]*dbasis[0*numnodes+j] 
								);
					/*term: u dvdx + v dudx + 2v dvdy -> u dvdx + 2v dvdy*/
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += D_scalar*thickness*basis[i]*(
								vx*basis[i]*dbasis[0*numnodes+j]
								+ 2*vy*basis[i]*dbasis[1*numnodes+j]
								);
				}
			}
		}else{
			_error_("Given isconvection = " << isconvection << " is not supported.");
		}

		/*Diffusion term: */
		factor = gauss->weight*Jdet*Ah*thickness*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[2*i*2*numnodes+2*j] += factor*(
							dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j]
							);
				Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*(
							dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j]
							);
			}
		}

		/*NOTE: For stability of stiffness matrix, the friction and coriolis terms are moved to Pvector section*/
		/*Friction term: Cd |u| (u, v)*/
		if(false){
			factor = gauss->weight*Jdet*Cd*vel*dt;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*term: Cd |u| u*/
					Ke->values[2*i*2*numnodes+2*j]       += factor*basis[i]*basis[j];
					/*term: Cd |u| v*/
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*basis[i]*basis[j];
				}
			}
		}

		/*Coriolis term: f (v, -u) */
		if(false){
			factor = gauss->weight*Jdet*dt*thickness*coriolis_freq;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*term: (f D v) */
					Ke->values[2*i*2*numnodes+2*j+1]   += -factor*basis[i]*basis[j];
					/*term: (-f D u)*/
					Ke->values[(2*i+1)*2*numnodes+2*j] += +factor*basis[i]*basis[j];
				}
			}
		}
		
		if(stabilization==1){ /* Artifical diffusivity {{{*/
			/*Artificial diffusivity: */
			if (true){
				D_scalar = gauss->weight*Jdet*dt*thickness;
				factor = D_scalar*h/2.0;

				/*Second-order tensor*/
				D[0][0]=factor*fabs(vx); D[0][1]=0.0;
				D[1][0]=0.0;             D[1][1]=factor*fabs(vy);
			}else{
				D_scalar = gauss->weight*Jdet*dt*thickness_avg;
				vel = sqrt(vx_avg*vx_avg + vy_avg*vy_avg);
				factor = D_scalar*h/2.0/vel;

				/*Second-order tensor*/
				D[0][0]=factor*fabs(vx_avg); D[0][1]=0.0;
				D[1][0]=0.0;                 D[1][1]=factor*fabs(vy_avg);
			}

			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[2*i*2*numnodes+2*j] += (
								dbasis[0*numnodes+i] *(D[0][0]*dbasis[0*numnodes+j] + D[0][1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1][0]*dbasis[0*numnodes+j] + D[1][1]*dbasis[1*numnodes+j]) 
								);

					Ke->values[(2*i+1)*2*numnodes+2*j+1] += (
								dbasis[0*numnodes+i] *(D[0][0]*dbasis[0*numnodes+j] + D[0][1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1][0]*dbasis[0*numnodes+j] + D[1][1]*dbasis[1*numnodes+j]) 
								);
				}
			}
		}/* }}} */
		else if(stabilization==2){ /*TODO: Stream-upwind scheme {{{ */
			/*Stream-upwind scheme: */
			/*TODO: Stream-upwind scheme should be fixed in future. This scheme is not working*/
			if(false){
				D_scalar = gauss->weight*Jdet*dt*thickness_avg;
				vel = sqrt(vx_avg*vx_avg* + vy_avg*vy_avg) + 1e-14;
				factor = D_scalar*h/(2.0*vel);

				/*Second-order tensor*/
				D[0][0]=factor*vx_avg*vx_avg;  D[0][1]=factor*vx_avg*vy_avg;
				D[1][0]=factor*vy_avg*vx_avg;  D[1][1]=factor*vy_avg*vy_avg;
			}else{
				D_scalar = gauss->weight*Jdet*dt*thickness;
				vel = sqrt(vx*vx* + vy*vy) + 1e-14;
				factor = D_scalar*h/(2.0*vel);

				/*Second-order tensor*/
				D[0][0]=factor*vx*vx;  D[0][1]=factor*vx*vy;
				D[1][0]=factor*vy*vx;  D[1][1]=factor*vy*vy;
			}

			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[2*i*2*numnodes+2*j] += (
								dbasis[0*numnodes+i] *(D[0][0]*dbasis[0*numnodes+j] + D[0][1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1][0]*dbasis[0*numnodes+j] + D[1][1]*dbasis[1*numnodes+j]) 
								);

					Ke->values[(2*i+1)*2*numnodes+2*j+1] += (
								dbasis[0*numnodes+i] *(D[0][0]*dbasis[0*numnodes+j] + D[0][1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1][0]*dbasis[0*numnodes+j] + D[1][1]*dbasis[1*numnodes+j]) 
								);
				}
			}
		}/* }}} */
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* BasalforcingsLaddieMomentumAnalysis::CreatePVectorCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int			stabilization,dim,domaintype;
	int         melt_style,point1;
	int         isentrainment;
	bool        mainlyfloating;
	IssmDouble  fraction1,fraction2;
	IssmDouble  Jdet,dt,intrusiondist;
	IssmDouble  vx,vy,vel,dvxdx,dvydy,xi,h,tau;
	IssmDouble  thickness; /*Plume thickness (or depth) [m]*/
	IssmDouble  dthk[2];
	IssmDouble  dthkdx, dthkdy;
	IssmDouble  dzb[2];
	IssmDouble  dzbdx, dzbdy;
	IssmDouble  g; /* gravitational acceleration [m s-2]*/
	IssmDouble  rho0; /*seawater density [kg m-3]*/
	IssmDouble  ga; /*ambient gravitational acceleration [m s-2*/
	IssmDouble  dentr; /*dentrainment rate*/
	IssmDouble  Cd;
	IssmDouble  f_cori;
	IssmDouble  ddrhoa[2];
	IssmDouble  ddrhoadx, ddrhoady;
	IssmDouble  gllevelset,phi=1.;
	IssmDouble *xyz_list = NULL;
	Gauss      *gauss     = NULL;

	IssmDouble factor_transient;
	IssmDouble factor_buoyancy;
	IssmDouble factor_pgradient;
	IssmDouble factor_dentr;
	IssmDouble factor_Cd, factor_cori;
	IssmDouble factor;

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
	ElementVector *pe    = element->NewElementVector();
	IssmDouble    *basis = xNew<IssmDouble>(numnodes);
	IssmDouble    *dbasis= xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationMomentumEnum);
	element->FindParam(&g,ConstantsGEnum);
	element->FindParam(&rho0,MaterialsRhoSeawaterEnum);
	element->FindParam(&Cd,BasalforcingsLaddieCdEnum);
	element->FindParam(&f_cori,BasalforcingsLaddieCoriolisFrequencyEnum);

	Input* thickness_input = element->GetInput(BasalforcingsLaddieThicknessOldEnum); _assert_(thickness_input);
	Input* zb_input        = element->GetInput(BaseEnum); _assert_(zb_input);
	Input* vx_input        = element->GetInput(BasalforcingsLaddieVxOldEnum); _assert_(vx_input);
	Input* vy_input        = element->GetInput(BasalforcingsLaddieVyOldEnum); _assert_(vy_input);
	Input* drhoa_input     = element->GetInput(BasalforcingsLaddieDRhoEnum); _assert_(drhoa_input);
	Input* ga_input        = element->GetInput(BasalforcingsLaddieAmbientGEnum); _assert_(ga_input);
	Input* dentr_input     = element->GetInput(BasalforcingsLaddieDEntrainmentRateEnum); _assert_(dentr_input);
	h=element->CharacteristicLength();

	/* Start looping on the number of gaussian points: */
	gauss=element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/*Get input values*/
		thickness_input->GetInputValue(&thickness,gauss);
		ga_input->GetInputValue(&ga,gauss);
		dentr_input->GetInputValue(&dentr,gauss);

		drhoa_input->GetInputDerivativeValue(&ddrhoa[0],xyz_list,gauss);
		thickness_input->GetInputDerivativeValue(&dthk[0],xyz_list,gauss);
		zb_input->GetInputDerivativeValue(&dzb[0],xyz_list,gauss);

		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=sqrt(vx*vx+vy*vy);

		ddrhoadx=ddrhoa[0];
		ddrhoady=ddrhoa[1];

		dthkdx=dthk[0];
		dthkdy=dthk[1];

		dzbdx=dzb[0];
		dzbdy=dzb[1];

		factor_transient= Jdet*gauss->weight;
		factor_buoyancy = Jdet*gauss->weight*dt*(-g*thickness*thickness/2/rho0);
		factor_pgradient= Jdet*gauss->weight*dt*(ga*thickness);
		//factor_dentr = 0.0;
		factor_dentr    = Jdet*gauss->weight*dt*dentr;
		factor_Cd       = Jdet*gauss->weight*dt*Cd*vel;
		factor_cori     = Jdet*gauss->weight*dt*thickness;

		for(int i=0;i<numnodes;i++){
			pe->values[i*2+0]+=basis[i]*(factor_transient*thickness*vx + factor_buoyancy*ddrhoadx + factor_pgradient*(dzbdx-dthkdx) - factor_dentr*vx + factor_Cd*vy - factor_Cd*vx);
			pe->values[i*2+1]+=basis[i]*(factor_transient*thickness*vy + factor_buoyancy*ddrhoady + factor_pgradient*(dzbdy-dthkdy) - factor_dentr*vy - factor_Cd*vx - factor_Cd*vy);
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return pe;
}/*}}}*/

/*Flux correction transport*/
//TODO: Flux corrected transport (FCT) is under developement.
ElementMatrix* BasalforcingsLaddieMomentumAnalysis::CreateFctKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         stabilization;
	int         domaintype,dim;
	IssmDouble  Jdet,D_scalar,h,factor;
	IssmDouble  dt; /*time-stepping*/
	IssmDouble  thickness;
	IssmDouble  vel,vx,vy,dvxdx,dvydy;
	IssmDouble  melt, entrain; /*unit m s-1*/
	IssmDouble  Cd;
	IssmDouble  dvx[2],dvy[2];
	IssmDouble  dthk[2];
	IssmDouble  dthkdx, dthkdy;
	IssmDouble  D[4];
	IssmDouble *xyz_list = NULL;
	IssmDouble  g; /*gravitational acceleration [m s-2]*/

	IssmDouble rho0; /*Sea water density [kg m-3]*/
	IssmDouble Kh; /*horizontal viscosity*/
	IssmDouble coriolis_freq; /*Coriolis frequency [s-1]*/

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
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationEnum);
	element->FindParam(&coriolis_freq,BasalforcingsLaddieCoriolisFrequencyEnum);
	element->FindParam(&Kh,BasalforcingsLaddieHorizontalViscosityEnum);
	element->FindParam(&Cd,BasalforcingsLaddieCdEnum);
	element->FindParam(&rho0,MaterialsRhoSeawaterEnum);
	element->FindParam(&g,ConstantsGEnum);
	element->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);

	/*Plume thickness D and depth averaged horizontal velocity*/
	Input* thickness_input = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thickness_input);
	Input* vx_input        = element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
	Input* vy_input        = element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);
	Input* rhoa_input      = element->GetInput(BasalforcingsLaddieDRhoEnum); _assert_(rhoa_input);
	Input* ga_input        = element->GetInput(BasalforcingsLaddieAmbientGEnum); _assert_(ga_input);

	h = element->CharacteristicLength();

	/* Start looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		
		/*Prepare inputs: */
		thickness_input->GetInputValue(&thickness,gauss);
		thickness_input->GetInputDerivativeValue(&dthk[0],xyz_list,gauss);
		dthkdx = dthk[0];
		dthkdy = dthk[1];

		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);

		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

		vel = sqrt(pow(vx,2.0) + pow(vy,2.0));;

		/*Convection term: */
		int isconvection=1;
		D_scalar=gauss->weight*Jdet*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				if(isconvection==0){
					/*term: 2u dudx + u dvdy + v dudy -> 2u dudx + v dudy*/
					Ke->values[2*i*2*numnodes+2*j] += -D_scalar*thickness*(
								+ 2*vx*basis[j]*dbasis[0*numnodes+i] 
								+ vy*basis[j]*dbasis[1*numnodes+i]
								);
					/*term: 2u dudx + u dvdy + v dudy -> u dvdy*/
					Ke->values[2*i*2*numnodes+2*j+1] += -D_scalar*thickness*(
								vx*basis[j]*dbasis[1*numnodes+i] 
								);

					/*term: u dvdx + v dudx + 2v dvdy -> v dudx */
					Ke->values[(2*i+1)*2*numnodes+2*j] += -D_scalar*thickness*(
								vy*basis[j]*dbasis[0*numnodes+i] 
								);
					/*term: u dvdx + v dudx + 2v dvdy -> u dvdx + 2v dvdy*/
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += D_scalar*thickness*(
								vx*basis[j]*dbasis[0*numnodes+i]
								+ 2*vy*basis[j]*dbasis[1*numnodes+i]
								);
				}else{
					/*term: 2u dudx + u dvdy + v dudy -> 2u dudx + v dudy*/
					Ke->values[2*i*2*numnodes+2*j] += -D_scalar*(
								+ 2*vx*basis[i]*(thickness*dbasis[0*numnodes+j] + dthkdx*basis[j])
								+ vy*basis[i]*(thickness*dbasis[1*numnodes+j] + dthkdy*basis[j])
								);
					/*term: 2u dudx + u dvdy + v dudy -> u dvdy*/
					Ke->values[2*i*2*numnodes+2*j+1] += -D_scalar*(
								basis[i]*vx*(thickness*dbasis[1*numnodes+j] + dthkdy*basis[j])
								);

					/*term: u dvdx + v dudx + 2v dvdy -> v dudx */
					Ke->values[(2*i+1)*2*numnodes+2*j] += -D_scalar*(
								basis[i]*vy*(thickness*dbasis[0*numnodes+j] + dthkdx*basis[j])
								);
					/*term: u dvdx + v dudx + 2v dvdy -> u dvdx + 2v dvdy*/
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += -D_scalar*(
								vx*basis[i]*(thickness*dbasis[0*numnodes+j] + dthkdx*basis[j])
								+ 2*basis[i]*vy*(thickness*dbasis[1*numnodes+j] + dthkdy*basis[j])
								);
				}
			}
		}

		/*Diffusion term: */
		factor = gauss->weight*Jdet*Kh*thickness*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[2*i*2*numnodes+2*j]       += -factor*(
							dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
							);
				Ke->values[(2*i+1)*2*numnodes+2*j+1] += -factor*(
							dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
							);
			}
		}

		/*Friction term: Cd |u| (u, v)*/
		factor = gauss->weight*Jdet*Cd*vel*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/*term: Cd |u| u*/
				Ke->values[2*i*2*numnodes+2*j]       += factor*basis[i]*basis[j];
				/*term: Cd |u| v*/
				Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*basis[i]*basis[j];
			}
		}

		/*Coriolis term: f (v, -u) */
		factor = gauss->weight*Jdet*coriolis_freq*thickness*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/*term: (f D v) */
				Ke->values[2*i*2*numnodes+2*j+1]   += -factor*basis[i]*basis[j];
				/*term: (-f D u)*/
				Ke->values[(2*i+1)*2*numnodes+2*j] += +factor*basis[i]*basis[j];
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementMatrix* BasalforcingsLaddieMomentumAnalysis::CreateMassMatrix(Element* element){/*{{{*/
	/* Check if ice in element */
	if(!element->IsOceanInElement() || !element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet;
	IssmDouble  thickness;
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Me     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input *thickness_input=element->GetInput(BasalforcingsLaddieThicknessEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		
		thickness_input->GetInputValue(&thickness,gauss);

		D=gauss->weight*Jdet*thickness;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Me->values[2*i*2*numnodes+2*j] += D*basis[i]*basis[j];
				Me->values[(2*i+1)*2*numnodes+2*j+1] += D*basis[i]*basis[j];
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Me;
}/*}}}*/
ElementVector* BasalforcingsLaddieMomentumAnalysis::CreateFctPVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement() || !element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int			stabilization,dim,domaintype;
	int         melt_style,point1;
	int         isentrainment;
	bool        mainlyfloating;
	IssmDouble  fraction1,fraction2;
	IssmDouble  Jdet,dt,intrusiondist;
	IssmDouble  vx,vy,vel,dvxdx,dvydy,xi,h,tau;
	IssmDouble  thickness; /*Plume thickness (or depth) [m]*/
	IssmDouble  dthk[2];
	IssmDouble  dthkdx, dthkdy;
	IssmDouble  dzb[2];
	IssmDouble  dzbdx, dzbdy;
	IssmDouble  g; /* gravitational acceleration [m s-2]*/
	IssmDouble  rho0; /*seawater density [kg m-3]*/
	IssmDouble  ga; /*ambient gravitational acceleration [m s-2*/
	IssmDouble  Cd;
	IssmDouble  ddrhoa[2];
	IssmDouble  ddrhoadx, ddrhoady;
	IssmDouble  gllevelset,phi=1.;
	IssmDouble *xyz_list = NULL;
	Gauss      *gauss     = NULL;

	IssmDouble factor_transient;
	IssmDouble factor_buoyancy;
	IssmDouble factor_pgradient;
	IssmDouble factor;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector *pe    = element->NewElementVector();
	IssmDouble    *basis = xNew<IssmDouble>(numnodes);
	IssmDouble    *dbasis= xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationMomentumEnum);

	Input* thickness_input= element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thickness_input);
	Input* ga_input= element->GetInput(BasalforcingsLaddieAmbientGEnum); _assert_(ga_input);
	Input* drhoa_input= element->GetInput(BasalforcingsLaddieDRhoEnum); _assert_(drhoa_input);

	Input* zb_input= element->GetInput(BasalforcingsLaddieAmbientGEnum); _assert_(zb_input);
	Input* vx_input   = element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
	Input* vy_input   = element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);
	Input* melt_input = element->GetInput(BasalforcingsLaddieMeltingRateEnum);  _assert_(melt_input);
	Input* entr_input = element->GetInput(BasalforcingsLaddieEntrainmentRateEnum);  _assert_(entr_input);

	/* Start  looping on the number of gaussian points: */
	gauss = element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/*Get input values*/
		thickness_input->GetInputValue(&thickness,gauss);
		ga_input->GetInputValue(&ga,gauss);

		drhoa_input->GetInputDerivativeValue(&ddrhoa[0],xyz_list,gauss);
		thickness_input->GetInputDerivativeValue(&dthk[0],xyz_list,gauss);
		zb_input->GetInputDerivativeValue(&dzb[0],xyz_list,gauss);

		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);

		ddrhoadx=ddrhoa[0];
		ddrhoady=ddrhoa[1];

		dthkdx=dthk[0];
		dthkdy=dthk[1];

		dzbdx=dzb[0];
		dzbdy=dzb[1];

		factor_transient= Jdet*gauss->weight;
		factor_buoyancy = Jdet*gauss->weight*dt*(-g*thickness*thickness/2/rho0);
		factor_pgradient= Jdet*gauss->weight*dt*(ga*thickness);

		for(int i=0;i<numnodes;i++){
			pe->values[i*2+0]+=(factor_transient*thickness*vx + factor_buoyancy*ddrhoadx + factor_pgradient*(dzbdx-dthkdx))*basis[i];
			pe->values[i*2+1]+=(factor_transient*thickness*vy + factor_buoyancy*ddrhoady + factor_pgradient*(dzbdy-dthkdy))*basis[i];
		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return pe;
}/*}}}*/
void           BasalforcingsLaddieMomentumAnalysis::FctKMatrix(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,FemModel* femmodel){/*{{{*/

	int domaintype;
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum) _error_("domain type "<<EnumToStringx(domaintype)<<" not supported yet");

	/*Output*/
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;

	/*Initialize Jacobian Matrix*/
	AllocateSystemMatricesx(&Kff,&Kfs,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(Object* & object : femmodel->elements->objects){
		Element*       element = xDynamicCast<Element*>(object);
		ElementMatrix* Ke     = this->CreateFctKMatrix(element);
		if(Ke) Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}
	Kff->Assemble();
	Kfs->Assemble();

	/*Assign output pointer*/
	*pKff=Kff;
	if(pKfs){
		*pKfs=Kfs;
	}
	else{
		delete Kfs;
	}
}/*}}}*/
void           BasalforcingsLaddieMomentumAnalysis::MassMatrix(Matrix<IssmDouble>** pMff,FemModel* femmodel){/*{{{*/

	/*Initialize Mass matrix*/
	Matrix<IssmDouble> *Mff = NULL;
	AllocateSystemMatricesx(&Mff,NULL,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(Object* & object : femmodel->elements->objects){
		Element*       element = xDynamicCast<Element*>(object);
		ElementMatrix* MLe     = this->CreateMassMatrix(element);
		if(MLe){
			MLe->AddToGlobal(Mff);
		}
		delete MLe;
	}
	Mff->Assemble();

	/*Assign output pointer*/
	*pMff=Mff;
}/*}}}*/
void           BasalforcingsLaddieMomentumAnalysis::FctPVector(Vector<IssmDouble>** ppf,FemModel* femmodel){/*{{{*/

	/*Output*/
	Vector<IssmDouble>* pf = NULL;

	/*Initialize P vector*/
	AllocateSystemMatricesx(NULL,NULL,NULL,&pf,femmodel);

	/*Create and assemble matrix*/
	for(Object* & object : femmodel->elements->objects){
		Element*       element = xDynamicCast<Element*>(object);
		ElementVector* pe      = this->CreateFctPVector(element);
		if(pe) pe->AddToGlobal(pf);
		delete pe;
	}
	pf->Assemble();

	/*Assign output pointer*/
	*ppf=pf;
}/*}}}*/
