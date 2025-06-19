#include "./BasalforcingsLaddieHeatAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/TransientInput.h"

#define FINITEELEMENT P1Enum

/*Model processing*/
void BasalforcingsLaddieHeatAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

}/*}}}*/
void BasalforcingsLaddieHeatAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
}/*}}}*/
void BasalforcingsLaddieHeatAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  stabilization;
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");

	/*Check in 3d*/
	if(stabilization==3 && iomodel->domaintype==Domain3DEnum) _error_("DG 3d not implemented yet");

	/*Create Nodes either DG or CG depending on stabilization*/
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,BasalforcingsLaddieHeatAnalysisEnum,P1Enum,isamr);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  BasalforcingsLaddieHeatAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void BasalforcingsLaddieHeatAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    stabilization;
	int    finiteelement;

	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");
	finiteelement=P1Enum;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}
}/*}}}*/
void BasalforcingsLaddieHeatAnalysis::UpdateParameters(Parameters* parBasalforcingsPlumeAnalysisEnumameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           BasalforcingsLaddieHeatAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BasalforcingsLaddieHeatAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BasalforcingsLaddieHeatAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BasalforcingsLaddieHeatAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BasalforcingsLaddieHeatAnalysis::CreateKMatrix(Element* element){/*{{{*/

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
ElementVector* BasalforcingsLaddieHeatAnalysis::CreatePVector(Element* element){/*{{{*/

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
void           BasalforcingsLaddieHeatAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,BasalforcingsLaddieTEnum);
}/*}}}*/
void           BasalforcingsLaddieHeatAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           BasalforcingsLaddieHeatAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Only update if on base*/
	if(!element->IsOnBase() || !element->IsIceInElement() || !element->IsOceanInElement()) return;

	int         i,dim,domaintype;
	int         numnodes;
	int        *doflist = NULL;
	IssmDouble *xyz_list;
	Element    *basalelement;

	element->FindParam(&domaintype,DomainTypeEnum);
	/*Get basal element*/
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim=2;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!element->IsOnBase()){xDelete<IssmDouble>(xyz_list); return;}
			basalelement=element->SpawnBasalElement();
			dim=2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	numnodes = basalelement->GetNumberOfNodes();
	int         numdof = numnodes*dim;
	IssmDouble *Tnew   = xNew<IssmDouble>(numnodes);

	/*Fetch dof list and allocate solution vector*/
	basalelement->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numnodes;i++){
		Tnew[i]=solution[doflist[i]];

		/*Check solution*/
		if(xIsNan<IssmDouble>(Tnew[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(Tnew[i])) _error_("Inf found in solution vector");
	}
	element->AddBasalInput(BasalforcingsLaddieTEnum,Tnew,element->GetElementType());

	xDelete<int>(doflist);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(Tnew);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           BasalforcingsLaddieHeatAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Deal with ocean constraint*/
	SetActiveNodesLSMx(femmodel);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);
		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);
		IssmDouble *Ta        = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&ls_active[0],IceMaskNodeActivationEnum);
		element->GetInputListOnNodes(&Ta[0],BasalforcingsLaddieAmbientTemperatureEnum);

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]<0. && ls_active[in]==1.){
				node->Activate();
			}
			else{
				/*Apply plume salt as zero value along grounding line*/
				node->Deactivate();
				if(true) node->ApplyConstraint(0,Ta[in]);
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(ls_active);
		xDelete<IssmDouble>(Ta);
	}
}/*}}}*/

/*Heat tranport analysis*/
ElementMatrix* BasalforcingsLaddieHeatAnalysis::CreateKMatrixCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement()) return NULL;

	/*Intermediaries */
	int         stabilization;
	int         domaintype,dim;
	IssmDouble  Jdet,D_scalar,dt,h,factor;
	IssmDouble  thickness; /*plume thickness*/
	IssmDouble  dthk[2];
	IssmDouble  dthkdx, dthkdy;
	IssmDouble  vel,vx,vy,dvxdx,dvydy;
	IssmDouble  thickness_avg, vx_avg, vy_avg;
	IssmDouble  xi,tau;
	IssmDouble  dvx[2],dvy[2];
	IssmDouble  D[2][2];
	IssmDouble *xyz_list = NULL;

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
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationEnum);

	element->FindParam(&Kh,BasalforcingsLaddieHorizontalDiffusivityEnum);

	/*Plume thickness D and depth averaged horizontal velocity*/
	Input* thickness_input = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thickness_input);
	Input* vx_input = element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
	Input* vy_input = element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);

	h = element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		/*Prepare inputs: */
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
		dvydy=dvy[1];

		/*Transient term*/
		/* d (DT) / dt = [D(t+1)T(t+1) - D(t)T(t)]/ dt*/
		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D_scalar*thickness*basis[i]*basis[j];
			}
		}

		/*Diffusion term: */
		D_scalar=gauss->weight*Jdet*dt*Kh*thickness;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/* thickness * dphi_i * dphi_j * T_j */
				Ke->values[i*numnodes+j] += D_scalar*(
							+ dbasis[0*numnodes+i]*dbasis[0*numnodes+j] 
							+ dbasis[1*numnodes+i]*dbasis[1*numnodes+j]
							);
			}
		}

		/*Advection term: */
		D_scalar=gauss->weight*Jdet*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				if(true){
					/*NOTE: d(u D T)/dx = u D dTdx + u T dDdx +  D T dudx =  u D dTdx + (u dDdx + D dudx) T*/
					
					/* (u dDdx + D dudx) T */
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*((vx*dthkdx + vy*dthkdy) + (thickness*dvxdx + thickness*dvydy));
					/* u D dTdx */
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*thickness*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
				}
				else{
					Ke->values[i*numnodes+j] += D_scalar*thickness*basis[i]*basis[j]*(dvxdx + dvydy);

					Ke->values[i*numnodes+j] += D_scalar*thickness*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
				}
			}
		}

		/*Prepare parameters for advection stabilization scheme*/
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) D[i][j]=0.0;

		if(stabilization==0){/*{{{*/
		}/*}}}*/
		else if(stabilization==1){/* Artificial diffusion {{{*/
			/*Artifical diffusion: */
			vx_input->GetInputValue(&vx,gauss);
			vy_input->GetInputValue(&vy,gauss);
			thickness_input->GetInputValue(&thickness,gauss);
			//vx_input->GetInputAverage(&vx);
			//vy_input->GetInputAverage(&vy);
			//thickness_input->GetInputAverage(&thickness);

			factor=D_scalar*h/2.0*thickness;
			D[0][0]=factor*fabs(vx); D[0][1]=0.0;
			D[1][0]=0.0;             D[1][1]=factor*fabs(vy);

			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += (
								dbasis[0*numnodes+i] *(D[0][0]*dbasis[0*numnodes+j] + D[0][1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1][0]*dbasis[0*numnodes+j] + D[1][1]*dbasis[1*numnodes+j])
								);
				}
			}
		}/*}}}*/
		else if(stabilization==2){/* Stream upwind {{{*/
			_assert_(dim==2);
			vx_input->GetInputAverage(&vx);
			vy_input->GetInputAverage(&vy);
			thickness_input->GetInputAverage(&thickness);

			/*Streamline upwind*/
			vel=sqrt(vx*vx+vy*vy)+1.e-14;
			tau=h/(2*vel);

			factor = dt*gauss->weight*Jdet*tau*thickness;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i])*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j]);
				}
			}
		}/*}}}*/
		else if(stabilization==5){/* Stream upwind Petrov Galerkin (SUPG) {{{*/
			/*SUPG*/
			if(dim!=2) _error_("Stabilization "<<stabilization<<" not supported yet for dim != 2");
			vx_input->GetInputAverage(&vx);
			vy_input->GetInputAverage(&vy);
			vel=sqrt(vx*vx+vy*vy)+1.e-8;
			//xi=0.3130;
			xi=1;
			tau=xi*h/(2*vel);
			//tau=dt/6; // as implemented in Ua

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
		else {/*{{{*/
			_error_("Stabilization "<<stabilization<<" not supported yet");
		}/*}}}*/
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* BasalforcingsLaddieHeatAnalysis::CreatePVectorCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int			stabilization,dim,domaintype;
	int         melt_style,point1;
	bool        mainlyfloating;
	IssmDouble  fraction1,fraction2;
	IssmDouble  Jdet,dt,intrusiondist;
	IssmDouble  thickness; /*plume thickness (or depth) [m]*/
	IssmDouble  T; /*plume temperature*/
	IssmDouble  Ta; /*ambient temperature*/
	IssmDouble  Tb; /*temprature at ice shelf base*/
	IssmDouble  gammaT; /*heat exchange coefficient*/
	IssmDouble  vx,vy,vel,dvxdx,dvydy,xi,h,tau;
	IssmDouble  melt_rate, entr_rate; /*melting rate and entrainment rate [unit: m s-1]*/
	IssmDouble  g; /* gravitational acceleration [m s-1]*/
	IssmDouble  dvx[2],dvy[2];
	IssmDouble  gllevelset,phi=1.;
	IssmDouble* xyz_list = NULL;
	Gauss*      gauss     = NULL;

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
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis= xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);
	element->FindParam(&stabilization,BasalforcingsLaddieStabilizationEnum);

	Input* thickness_input= element->GetInput(BasalforcingsLaddieThicknessOldEnum);   _assert_(thickness_input);
	Input* vx_input       = element->GetInput(BasalforcingsLaddieVxEnum);		  _assert_(vx_input);
	Input* vy_input       = element->GetInput(BasalforcingsLaddieVyEnum);	     _assert_(vy_input);
	Input* T_input        = element->GetInput(BasalforcingsLaddieTEnum);   _assert_(T_input);
	Input* Tb_input       = element->GetInput(BasalforcingsLaddieTbEnum);   _assert_(Tb_input);
	Input* melt_input     = element->GetInput(BasalforcingsLaddieMeltingRateEnum); _assert_(melt_input);
	Input* entr_input     = element->GetInput(BasalforcingsLaddieEntrainmentRateEnum); _assert_(melt_input);
	Input* gammaT_input   = element->GetInput(BasalforcingsLaddieGammaTEnum); _assert_(gammaT_input);
	/*Ambient temperature*/
	Input* Ta_input       = element->GetInput(BasalforcingsLaddieAmbientTemperatureEnum); _assert_(Ta_input);
	h=element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	gauss = element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/*Get input values*/
		thickness_input->GetInputValue(&thickness,gauss); /*plumem thickness*/
		T_input->GetInputValue(&T,gauss); /*plume temperature*/
		Ta_input->GetInputValue(&Ta,gauss); /*ambient temperature*/
		Tb_input->GetInputValue(&Tb,gauss); /*temperature at ice shelf base*/
		melt_input->GetInputValue(&melt_rate,gauss);
		entr_input->GetInputValue(&entr_rate,gauss);
		gammaT_input->GetInputValue(&gammaT,gauss);

		IssmDouble factor = Jdet*gauss->weight*(thickness*T+dt*(entr_rate*Ta + melt_rate*Tb - gammaT*(T-Tb)));
		for(int i=0;i<numnodes;i++){
			pe->values[i]+=factor*basis[i];
		}

		if(stabilization==5){ //SUPG {{{
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
			factor = Jdet*gauss->weight*(thickness*T+dt*(entr_rate*gammaT + melt_rate*Tb - gammaT*(T-Tb)));
			for(int i=0;i<numnodes;i++){
				pe->values[i]+=factor*(tau*vx*dbasis[0*numnodes+i]+tau*vy*dbasis[1*numnodes+i]);
			}
			/*Force vector - part 3*/
			for(int i=0;i<numnodes;i++){
				pe->values[i]+=factor*(tau*basis[i]*dvxdx+tau*basis[i]*dvydy);
			}
		} //}}}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return pe;
}/*}}}*/
