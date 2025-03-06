#include "./BasalforcingsPlumeAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/TransientInput.h"

#define FINITEELEMENT P1Enum

/*Model processing*/
void BasalforcingsPlumeAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int stabilization;
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");

	/*Do not add constraints in DG,  they are weakly imposed*/
	if(stabilization!=3){
		IoModelToConstraintsx(constraints,iomodel,"md.basalforcings.spcthickness",MasstransportAnalysisEnum,FINITEELEMENT);
	}
}/*}}}*/
void BasalforcingsPlumeAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

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
void BasalforcingsPlumeAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  stabilization;
	iomodel->FindConstant(&stabilization,"md.masstransport.stabilization");

	/*Check in 3d*/
	if(stabilization==3 && iomodel->domaintype==Domain3DEnum) _error_("DG 3d not implemented yet");

	/*Create Nodes either DG or CG depending on stabilization*/
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	if(stabilization!=3){
		::CreateNodes(nodes,iomodel,BasalforcingsPlumeAnalysisEnum,FINITEELEMENT,isamr);
	}
	else{
		::CreateNodes(nodes,iomodel,BasalforcingsPlumeAnalysisEnum,P1DGEnum,isamr);
	}
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  BasalforcingsPlumeAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void BasalforcingsPlumeAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    stabilization,finiteelement;
	int    grdmodel;

	/*Fetch data needed: */
	iomodel->FindConstant(&stabilization,"md.basalforcings.stabilization");

	/*Finite element type*/
	finiteelement = FINITEELEMENT;
	if(stabilization==3){
		finiteelement = P1DGEnum;
		//finiteelement = P0DGEnum;
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

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.D",BasalforcingsPlumeDepthEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.vx",BasalforcingsPlumeVyEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.vy",BasalforcingsPlumeVxEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.temperature",BasalforcingsPlumeTEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.salinity",BasalforcingsPlumeSEnum,0);

	if(isgroundingline) 	iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
	/*Initialize ThicknessResidual input*/
	InputUpdateFromConstantx(inputs,elements,0.,ThicknessResidualEnum);

	/*Get what we need for ocean-induced basal melting*/
	bool isstochastic;
	int basalforcing_model;
	int melt_parameterization;
	iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
	iomodel->FindConstant(&isstochastic,"md.stochasticforcing.isstochasticforcing");
	iomodel->FindConstant(&melt_parameterization,"md.frontalforcings.parameterization");
	switch(basalforcing_model){
		case FloatingMeltRateEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.floatingice_melting_rate",BasalforcingsFloatingiceMeltingRateEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.perturbation_melting_rate",BasalforcingsPerturbationMeltingRateEnum,0.);
			if(isstochastic){
            iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
            iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.floatingice_melting_rate",BaselineBasalforcingsFloatingiceMeltingRateEnum);
         }
			break;
		case LinearFloatingMeltRateEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.perturbation_melting_rate",BasalforcingsPerturbationMeltingRateEnum,0.);
			break;
		case MismipFloatingMeltRateEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.meltrate_factor",BasalforcingsMeltrateFactorEnum);
			break;
		case MantlePlumeGeothermalFluxEnum:
			break;
		case SpatialLinearFloatingMeltRateEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_melting_rate",BasalforcingsSpatialDeepwaterMeltingRateEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_elevation",BasalforcingsSpatialDeepwaterElevationEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.upperwater_melting_rate",BasalforcingsSpatialUpperwaterMeltingRateEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.upperwater_elevation",BasalforcingsSpatialUpperwaterElevationEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.perturbation_melting_rate",BasalforcingsPerturbationMeltingRateEnum,0.);
			if(isstochastic){
            iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
            iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_melting_rate",BaselineBasalforcingsSpatialDeepwaterMeltingRateEnum);
         }
			break;
		case BasalforcingsPicoEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basin_id",BasalforcingsPicoBasinIdEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.overturning_coeff",BasalforcingsPicoOverturningCoeffEnum);
			break;
		case BasalforcingsIsmip6Enum:{
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basin_id",BasalforcingsIsmip6BasinIdEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.melt_anomaly",BasalforcingsIsmip6MeltAnomalyEnum,0.);

			/*Deal with tf...*/
			IssmDouble* array2d = NULL; int M,N,K; IssmDouble* temp = NULL;
			iomodel->FetchData(&temp,&M,&K,"md.basalforcings.tf_depths"); xDelete<IssmDouble>(temp);
			_assert_(M==1); _assert_(K>=1);
			for(int kk=0;kk<K;kk++){

				/*Fetch TF for this depth*/
				iomodel->FetchData(&array2d, &M, &N, kk, "md.basalforcings.tf");
				if(!array2d) _error_("md.basalforcings.tf not found in binary file");
				for(Object* & object : elements->objects){
					Element*  element = xDynamicCast<Element*>(object);
					if(iomodel->domaintype!=Domain2DhorizontalEnum && !element->IsOnBase()) continue;
					element->DatasetInputAdd(BasalforcingsIsmip6TfEnum,array2d,inputs,iomodel,M,N,1,BasalforcingsIsmip6TfEnum,kk);
				}
			xDelete<IssmDouble>(array2d);
			}
											  }
			break;
		case BeckmannGoosseFloatingMeltRateEnum:
			bool isthermalforcing;
			iomodel->FindConstant(&isthermalforcing,"md.basalforcings.isthermalforcing");
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.meltrate_factor",BasalforcingsMeltrateFactorEnum);
			if(isthermalforcing==0){
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_salinity",BasalforcingsOceanSalinityEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_temp",BasalforcingsOceanTempEnum);
			}
			else if(melt_parameterization!=FrontalForcingsRignotarmaEnum){
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_thermalforcing",ThermalForcingEnum);
			}
			break;
		case LinearFloatingMeltRatearmaEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basin_id",BasalforcingsLinearBasinIdEnum);
			if(isstochastic) iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
			break;
		default:
			_error_("Basal forcing model "<<EnumToStringx(basalforcing_model)<<" not supported yet");
	}

	if(!issmb){
		iomodel->FetchDataToInput(inputs,elements,"md.smb.mass_balance",SmbMassBalanceEnum);
	}
	if(stabilization==3){
		iomodel->FetchDataToInput(inputs,elements,"md.masstransport.spcthickness",MasstransportSpcthicknessEnum); //for DG, we need the spc in the element
	}
	if(stabilization==4){
		iomodel->FetchDataToInput(inputs,elements,"md.masstransport.spcthickness",MasstransportSpcthicknessEnum); //for FCT, we need the spc in the element (penlaties)
	}

	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}

	/*Initialize sea level cumulated sea level loads :*/
	iomodel->ConstantToInput(inputs,elements,0,AccumulatedDeltaIceThicknessEnum,P1Enum);
	iomodel->ConstantToInput(inputs,elements,0,OldAccumulatedDeltaIceThicknessEnum,P1Enum);

	/*for Ivins deformation model, initialize history of ice thickness changes:*/
	iomodel->FindConstant(&grdmodel,"md.solidearth.settings.grdmodel");
	if(grdmodel==IvinsEnum) inputs->SetTransientInput(TransientAccumulatedDeltaIceThicknessEnum,NULL,0);

}/*}}}*/
void BasalforcingsPlumeAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isFS",FlowequationIsFSEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.isfreesurface",MasstransportIsfreesurfaceEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.hydrostatic_adjustment",MasstransportHydrostaticAdjustmentEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.stabilization",MasstransportStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.min_thickness",MasstransportMinThicknessEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.penalty_factor",MasstransportPenaltyFactorEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.groundingline.intrusion_distance",GroundinglineIntrusionDistanceEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.masstransport.requested_outputs");
	parameters->AddObject(new IntParam(MasstransportNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(MasstransportRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.masstransport.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           BasalforcingsPlumeAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BasalforcingsPlumeAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BasalforcingsPlumeAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BasalforcingsPlumeAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BasalforcingsPlumeAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementMatrix* Ke = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			Ke = CreateKMatrixCG(basalelement);
			break;
		case P0DGEnum:
		case P1DGEnum:
			Ke = CreateKMatrixDG(basalelement);
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementMatrix* BasalforcingsPlumeAnalysis::CreateKMatrixCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement()) return NULL;

	/*Intermediaries */
	int        stabilization;
	int        domaintype,dim;
	IssmDouble Jdet,D_scalar,dt,h,factor;
	IssmDouble vel,vx,vy,dvxdx,dvydy;
	IssmDouble xi,tau;
	IssmDouble dvx[2],dvy[2];
	IssmDouble D[4];
	IssmDouble* xyz_list = NULL;

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
	element->FindParam(&stabilization,MasstransportStabilizationEnum);
	Input* vxaverage_input=element->GetInput(BasalforcingsPlumeVxEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=element->GetInput(BasalforcingsPlumeVyEnum); _assert_(vyaverage_input);
	/*Plume thickness D*/
	Input* thickness_input=element->GetInput(BasalforcingsPlumeDepthEnum); _assert_(thickness_input);

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
		vxaverage_input->GetInputValue(&vx,gauss);
		vxaverage_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		D_scalar=dt*gauss->weight*Jdet;
		if(dim==2){
			vyaverage_input->GetInputValue(&vy,gauss);
			vyaverage_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
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
		}
		else{
			dvxdx=dvx[0];
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += D_scalar*dvxdx*basis[i]*basis[j];
					Ke->values[i*numnodes+j] += D_scalar*vx*dbasis[0*numnodes+j]*basis[i];
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
				vxaverage_input->GetInputAverage(&vx);
				if(dim==2) vyaverage_input->GetInputAverage(&vy);
				D[0*dim+0]=h/2.0*fabs(vx);
				if(dim==2) D[1*dim+1]=h/2.0*fabs(vy);
				break;
			case 2:
				/*Streamline upwinding*/
				vxaverage_input->GetInputAverage(&vx);
				if(dim==1){
					vel=fabs(vx)+1.e-8;
				}
				else{
					vyaverage_input->GetInputAverage(&vy);
					vel=sqrt(vx*vx+vy*vy)+1.e-8;
				}
				tau=h/(2*vel);
				break;
			case 5:
				/*SUPG*/
				if(dim!=2) _error_("Stabilization "<<stabilization<<" not supported yet for dim != 2");
				vxaverage_input->GetInputAverage(&vx);
				vyaverage_input->GetInputAverage(&vy);
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				//xi=0.3130;
				xi=1;
				tau=xi*h/(2*vel);
				//tau=dt/6; // as implemented in Ua
				break;
			default:
				_error_("Stabilization "<<stabilization<<" not supported yet");
		}
		if(stabilization==1){
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
		}
		if(stabilization==2){
			/*Streamline upwind*/
			_assert_(dim==2);
			factor = dt*gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i])*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j]);
				}
			}
		}
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
ElementVector* BasalforcingsPlumeAnalysis::CreatePVector(Element* element){/*{{{*/

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
			pe = CreatePVectorDG(basalelement);
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
ElementVector* BasalforcingsPlumeAnalysis::CreatePVectorCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int			stabilization,dim,domaintype;
	int         melt_style,point1;
	bool        mainlyfloating;
	IssmDouble  fraction1,fraction2;
	IssmDouble  Jdet,dt,intrusiondist;
	IssmDouble  ms,mb,gmb,fmb,thickness,fmb_pert,gldistance;
	IssmDouble  vx,vy,vel,dvxdx,dvydy,xi,h,tau;
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
	element->FindParam(&melt_style,GroundinglineMeltInterpolationEnum);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,MasstransportStabilizationEnum);
	element->FindParam(&intrusiondist,GroundinglineIntrusionDistanceEnum);

	Input* gmb_input        = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);  _assert_(gmb_input);
	Input* fmb_input        = element->GetInput(BasalforcingsFloatingiceMeltingRateEnum);  _assert_(fmb_input);
	//Input* fmb_pert_input   = element->GetInput(BasalforcingsPerturbationMeltingRateEnum); _assert_(fmb_pert_input);
	Input* gllevelset_input = element->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
	Input* base_input = element->GetInput(BaseEnum);                    _assert_(base_input);
	Input* vx_input = element->GetInput(BasalforcingsPlumeVxEnum);		     _assert_(vxaverage_input);
	Input* vy_input = element->GetInput(BasalforcingsPlumeVyEnum);	        _assert_(vyaverage_input);
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

		ms_input->GetInputValue(&ms,gauss);
		gmb_input->GetInputValue(&gmb,gauss);
		fmb_input->GetInputValue(&fmb,gauss);
		//fmb_pert_input->GetInputValue(&fmb_pert,gauss);
		gllevelset_input->GetInputValue(&gllevelset,gauss);
		thickness_input->GetInputValue(&thickness,gauss);

		if(melt_style==SubelementMelt1Enum){
			if (phi>0.999999999) mb=gmb;
			else mb=(1-phi)*fmb+phi*gmb; // phi is the fraction of grounded ice so (1-phi) is floating
		}
		else if(melt_style==SubelementMelt2Enum){
			if(gllevelset>0.) mb=gmb;
			else mb=fmb;
		}
		else if(melt_style==NoMeltOnPartiallyFloatingEnum){
			if (phi<0.00000001){
				mb=fmb;//+fmb_pert;
			}
			else mb=gmb;
		}
		else if(melt_style==FullMeltOnPartiallyFloatingEnum){
			if (phi<0.99999999) mb=fmb;
			else mb=gmb;
		}
		else if(melt_style==IntrusionMeltEnum){
			Input* gldistance_input = element->GetInput(DistanceToGroundinglineEnum); _assert_(gldistance_input); 
			gldistance_input->GetInputValue(&gldistance,gauss);
			if(intrusiondist==0){
				if(gllevelset>0.) mb=gmb;
				else mb=fmb;
			}
			else if(gldistance>intrusiondist) {
				mb=gmb;
			}
			else if(gldistance<=intrusiondist && gldistance>0) {
				mb=fmb*(1-gldistance/intrusiondist); 
			}
			else{
				mb=fmb;
			}
		}
		else{
			_error_("melt interpolation "<<EnumToStringx(melt_style)<<" not implemented yet");
		}

		IssmDouble factor = Jdet*gauss->weight*(thickness+dt*(ms-mb));
		for(int i=0;i<numnodes;i++) pe->values[i]+=factor*basis[i];

		if(stabilization==5){ //SUPG
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			vxaverage_input->GetInputAverage(&vx);
			vyaverage_input->GetInputAverage(&vy);
			vxaverage_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vyaverage_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			vel=sqrt(vx*vx+vy*vy)+1.e-8;
			dvxdx=dvx[0];
			dvydy=dvy[1];
			//xi=0.3130;
			xi=1;
			tau=xi*h/(2*vel);
			//tau=dt/6; // as implemented in Ua

			/*Force vector - part 2*/
			factor = Jdet*gauss->weight*(thickness+dt*(ms-mb));
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

ElementMatrix* BasalforcingsPlumeAnalysis::CreateKMatrixMassbalance(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementMatrix* Ke = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			Ke = CreateKMatrixCGMassbalance(basalelement);
			break;
		case P0DGEnum:
		case P1DGEnum:
			Ke = CreateKMatrixDGMassbalance(basalelement);
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementMatrix* BasalforcingsPlumeAnalysis::CreateKMatrixCGMassbalance(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsOceanInElement()) return NULL;

	/*Intermediaries */
	int        stabilization;
	int        domaintype,dim;
	::IssmDouble Jdet,D_scalar,dt,h,factor;
	IssmDouble vel,vx,vy,dvxdx,dvydy;
	IssmDouble xi,tau;
	IssmDouble dvx[2],dvy[2];
	IssmDouble D[4];
	IssmDouble* xyz_list = NULL;

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
	element->FindParam(&stabilization,MasstransportStabilizationEnum);
	Input* vxaverage_input=element->GetInput(BasalforcingsPlumeVxEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=element->GetInput(BasalforcingsPlumeVyEnum); _assert_(vyaverage_input);
	/*Plume thickness D*/
	Input* thickness_input=element->GetInput(BasalforcingsPlumeDepthEnum); _assert_(thickness_input);

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
		vxaverage_input->GetInputValue(&vx,gauss);
		vxaverage_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		D_scalar=dt*gauss->weight*Jdet;
		if(dim==2){
			vyaverage_input->GetInputValue(&vy,gauss);
			vyaverage_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
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
		}
		else{
			dvxdx=dvx[0];
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += D_scalar*dvxdx*basis[i]*basis[j];
					Ke->values[i*numnodes+j] += D_scalar*vx*dbasis[0*numnodes+j]*basis[i];
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
				vxaverage_input->GetInputAverage(&vx);
				if(dim==2) vyaverage_input->GetInputAverage(&vy);
				D[0*dim+0]=h/2.0*fabs(vx);
				if(dim==2) D[1*dim+1]=h/2.0*fabs(vy);
				break;
			case 2:
				/*Streamline upwinding*/
				vxaverage_input->GetInputAverage(&vx);
				if(dim==1){
					vel=fabs(vx)+1.e-8;
				}
				else{
					vyaverage_input->GetInputAverage(&vy);
					vel=sqrt(vx*vx+vy*vy)+1.e-8;
				}
				tau=h/(2*vel);
				break;
			case 5:
				/*SUPG*/
				if(dim!=2) _error_("Stabilization "<<stabilization<<" not supported yet for dim != 2");
				vxaverage_input->GetInputAverage(&vx);
				vyaverage_input->GetInputAverage(&vy);
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				//xi=0.3130;
				xi=1;
				tau=xi*h/(2*vel);
				//tau=dt/6; // as implemented in Ua
				break;
			default:
				_error_("Stabilization "<<stabilization<<" not supported yet");
		}
		if(stabilization==1){
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
		}
		if(stabilization==2){
			/*Streamline upwind*/
			_assert_(dim==2);
			factor = dt*gauss->weight*Jdet*tau;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=factor*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i])*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j]);
				}
			}
		}
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

IssmDouble GetBasalforcingsFrictionVelocity(Element *element){/*{{{*/
	/*
	 Calculate the friction velocity Ustar defined in Jenkins et al. (2010)

	 Ustar = sqrt(Cd_top (vx^2 + vy^2 + Utide^2))

	 See also
	 Eq 13. in Lambert et al. (2023TC)

	 Output
	 * Ustar - element friction velocity
	 */

	/*Hard coding for specific parameters*/
	IssmDouble Cd_top=1.1e-3; // top drag coefficient 
	IssmDouble Utide=0.01; // tide velocity [m s-1]
	
	/*Input values*/
	IssmDouble vx, vy;

	/*Output values*/
	IssmDouble Ustar; // Friction velocity

	/* Get parameters */
	element->FindParam(&Cd_top,BasalforcingsPlumeCdTop);

	/* Get inputs */
	Input* vx_input =element->GetInput(BasalforcingsPlumeVxEnum); _assert_(vx_input);
	Input* vy_input =element->GetInput(BasalforcingsPlumeVyEnum); _assert_(vy_input);

	/*Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(1);
	/*Retrieve input value */
	gauss->GaussPoint(0); /* Maybe element value */

	vx_input->GetInputValue(&vx,gauss);
	vy_input->GetInputValue(&vy,gauss);

	/*Calculate fricition velocity*/
	Ustar = (Cd_top * (vx^2 + vy^2 + Utide^2))**0.5;

	return Ustar
}/*}}}*/
void GetHeatExchangeCoefficient(Element *element, IssmDouble *pgammaT, IssmDouble *pgammaS){/*{{{*/
	/*
	 Explain
	 Calculate elemental heat exchange coefficient for temperature (gamma_T) and salinity (gamma_S).

	 Output
	 * pgammaT
	 * pgammaS

	 See also Eq. 11 and Eq. 12 in Lambert et al. (2023TC).
	 */

	/*Initialize variable*/
	IssmDouble Ustar; // Friction velocity 
	IssmDouble Pr=13.8; /*Prandtl number*/
	IssmDouble Sc=2432; /*Schmidt number*/
	
	/*Get friction velocity*/
	Ustar = GetBasalforcingsFrictionVelocity(element);

	/*Calculate exchange coefficient for temperature*/
	pgammaT = Ustar/(2.12*log10(Ustar/D)
}/*}}}*/
