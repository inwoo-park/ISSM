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
		::CreateNodes(nodes,iomodel,BasalforcingsPlumeAnalysisEnum,FINITEELEMENT,isamr);
	}
	else{
		::CreateNodes(nodes,iomodel,BasalforcingsPlumeAnalysisEnum,P1DGEnum,isamr);
	}
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  BasalforcingsLaddieSaltAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void BasalforcingsLaddieSaltAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

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
void BasalforcingsLaddieSaltAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	/*Basalforcing plume model specials*/
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Kh",BasalforcingsLaddieHorizontalDiffusivityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Ah",BasalforcingsLaddieHorizontalViscosityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Dmin",BasalforcingsLaddieDepthMinEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Utide",BasalforcingsLddieTideVelocityEnum));

	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.isentrainment",BasalforcingsLaddieIsEntrainmentEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.Cd_mon",BasalforcingsLaddieCdMomentumEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.f_cori",BasalforcingsLaddieCoriolisFrequencyEnum));

	//iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.basal.requested_outputs");
	//parameters->AddObject(new IntParam(MasstransportNumRequestedOutputsEnum,numoutputs));
	//if(numoutputs)parameters->AddObject(new StringArrayParam(MasstransportRequestedOutputsEnum,requestedoutputs,numoutputs));
	//iomodel->DeleteData(&requestedoutputs,numoutputs,"md.masstransport.requested_outputs");

}/*}}}*/
