#include "./DamageEvolutionAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void DamageEvolutionAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters*/
	int damage_model;
	iomodel->FindConstant(&damage_model,"md.damage.model");
	if(damage_model != DamageEvolutionEnum) return;

	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.damage.elementinterp");

	/*Fetch parameters: */
	int stabilization;
	iomodel->FindConstant(&stabilization,"md.damage.stabilization");

	/*Do not add constraints in DG,  they are weakly imposed*/
	if(stabilization!=3){
		IoModelToConstraintsx(constraints,iomodel,"md.damage.spcdamage",DamageEvolutionAnalysisEnum,finiteelement);
	}

	/*FCT, constraints are imposed using penalties*/
	if(stabilization==4){
		constraints->ActivatePenaltyMethod(DamageEvolutionAnalysisEnum);
	}
}/*}}}*/
void DamageEvolutionAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Nothing for now*/

}/*}}}*/
void DamageEvolutionAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int finiteelement;

	iomodel->FindConstant(&finiteelement,"md.damage.elementinterp");
	::CreateNodes(nodes,iomodel,DamageEvolutionAnalysisEnum,finiteelement);
}/*}}}*/
int  DamageEvolutionAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void DamageEvolutionAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int finiteelement;
	bool   ismovingfront;

	iomodel->FindConstant(&finiteelement,"md.damage.elementinterp");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");

	/*Update elements: */
	iomodel->FetchData(1,"md.flowequation.element_equation");
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);

			/*Need to know the type of approximation for this element*/
			if(iomodel->Data("md.flowequation.element_equation")){
				inputs->SetInput(ApproximationEnum,counter,IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[i])));
			}
			counter++;
		}
	}
	iomodel->DeleteData(1,"md.flowequation.element_equation");

   /*First, reset all F to 0 */
	for(Object* & object : elements->objects){
		Element*       element = xDynamicCast<Element*>(object);
		element->SetElementInput(inputs,DamageFEnum,0.);
	}

	/*What input do I need to run my damage evolution model?*/
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VxMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VyMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VzMeshEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.damage.D",DamageDEnum);
	if (iomodel->domaintype==Domain3DEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.damage.D",DamageDbarEnum,0.0);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);

}/*}}}*/
void DamageEvolutionAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*Intermediaries*/
	int     damage_model;
	int     numoutputs;
	char**  requestedoutputs = NULL;
	iomodel->FindConstant(&damage_model,"md.damage.model");

	/*Now, do we really want DamageEvolution?*/
	if(damage_model!=DamageEvolutionEnum) return;

	/*retrieve some parameters: */
	parameters->AddObject(new IntParam(DamageModelEnum,damage_model));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.law",DamageLawEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.stabilization",DamageStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.max_damage",DamageMaxDamageEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.critical_damage",DamageCriticalDamageEnum));

	/*Requested outputs*/
	//NOTE: this section is duplicted at "modules/ModelProcessorx/CreateParameters.cpp".
	//iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.damage.requested_outputs");
	//parameters->AddObject(new IntParam(DamageEvolutionNumRequestedOutputsEnum,numoutputs));
	//if(numoutputs)parameters->AddObject(new StringArrayParam(DamageEvolutionRequestedOutputsEnum,requestedoutputs,numoutputs));
	//iomodel->DeleteData(&requestedoutputs,numoutputs,"md.damage.requested_outputs");

	/*Retrieve law dependent parameters: */
	int law;
	iomodel->FindConstant(&law,"md.damage.law");
	if (law==0){
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.stress_threshold",DamageStressThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.kappa",DamageKappaEnum));
	}
	else if (law>0){
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c1",DamageC1Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c2",DamageC2Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c3",DamageC3Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c4",DamageC4Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.stress_threshold",DamageStressThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.stress_ubound",DamageStressUBoundEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.kappa",DamageKappaEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.healing",DamageHealingEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.equiv_stress",DamageEquivStressEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.isdamage_exponent",DamageIsDamageExponentEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.isPeff",DamageIsPeffEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.iseffective_stress",DamageIsEffectiveStressEnum));
	}

	/* Parameters relevant for stress criterion */
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.equiv_stress_alpha",DamageEquivStressAlphaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.equiv_stress_beta",DamageEquivStressBetaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.equiv_stress_mu",DamageEquivStressMuEnum));

	/* Relevant for flow equation */
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isSIA",FlowequationIsSIAEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isSSA",FlowequationIsSSAEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isL1L2",FlowequationIsL1L2Enum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isMOLHO",FlowequationIsMOLHOEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isHO",FlowequationIsHOEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isFS",FlowequationIsFSEnum));

}/*}}}*/

/*Finite Element Analysis*/
void           DamageEvolutionAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInput(Element* element){/*{{{*/

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*Calculate damage evolution source term: */
	for (int i=0;i<numnodes;i++){

		/* healing could be handled here */

		/* no source term; damage handled in stress balance */
		f[i]=0.;
	}

	/*Add input*/
	element->AddInput(DamageFEnum,f,element->GetElementType());

	/*Clean up and return*/
	xDelete<IssmDouble>(f);
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInputArctan(Element* element){/*{{{*/
	IssmDouble c1, c2, stress_threshold, stress_ubound;
	IssmDouble damage;
	IssmDouble yts;
	IssmDouble principalDevStress1, principalDevStress2;
	IssmDouble tensileStress, compressiveStress;

	int domaintype, dim;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*retrieve parameters:*/
	element->FindParam(&c1,DamageC1Enum);
	element->FindParam(&c2,DamageC2Enum);
	element->FindParam(&yts,ConstantsYtsEnum);
	element->FindParam(&stress_threshold,DamageStressThresholdEnum);
	element->FindParam(&stress_ubound,DamageStressUBoundEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Get problem dimension*/
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented");
	}
	/*Compute stress tensor and Stress Max Principal: */
	element->ComputeDeviatoricStressTensor();

	Input* principalDevStress1_input = element->GetInput(DeviatoricStress1Enum);     _assert_(principalDevStress1_input);
	Input* principalDevStress2_input = element->GetInput(DeviatoricStress2Enum);     _assert_(principalDevStress2_input);

	Input* damage_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
	}

	/*Calculate damage evolution source term */
	Gauss* gauss=element->NewGauss();

	/* To keep arctan output (bounded by -pi/2 and pi/2) within the specified boundaries */
	c1 /= (PI/2);
	c2 /= (PI/2);
	/* To have per second output with per annum parameters */
	c1 /= yts;
	c2 /= yts;

	for (int i=0;i<numnodes;i++){
		f[i] = 0;

		gauss->GaussNode(element->GetElementType(),i);

		damage_input->GetInputValue(&damage,gauss);
		principalDevStress1_input->GetInputValue(&principalDevStress1,gauss);
		principalDevStress2_input->GetInputValue(&principalDevStress2,gauss);

		tensileStress     = sqrt(1.5*(pow(max(principalDevStress1, 0.), 2) + pow(max(principalDevStress2, 0.), 2)));
		compressiveStress = sqrt(1.5*(pow(min(principalDevStress1, 0.), 2) + pow(min(principalDevStress2, 0.), 2)));

		/* Calculate principal effective stresses */
		if(dim==2){
			f[i] = 0;
			if(tensileStress > stress_threshold)
				f[i] += c1*atan((tensileStress/stress_threshold - 1)/(1-damage));

			if(compressiveStress < stress_ubound)
				f[i] += c2*atan((compressiveStress/stress_ubound - 1)/(1-damage));
		}
		else{
			_error_("Only 2D is implemented.");
		}
	}

	/*Add input*/
	element->AddInput(DamageFEnum,f,element->GetElementType());

	/*Clean up and return*/
	xDelete<IssmDouble>(f);
	delete gauss;
}/*}}}*/

void           DamageEvolutionAnalysis::CreateDamageFInputExp(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble epsf,stress_threshold,eps0;
	IssmDouble damage,B,n,epseff;
	IssmDouble eps_xx,eps_yy,eps_xy,eps1,eps2,epstmp;
	int domaintype;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*retrieve parameters:*/
	element->FindParam(&epsf,DamageC1Enum);
	element->FindParam(&stress_threshold,DamageStressThresholdEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Compute stress tensor: */
	element->ComputeStrainRate();

	/*retrieve what we need: */
	Input* eps_xx_input  = element->GetInput(StrainRatexxEnum);     _assert_(eps_xx_input);
	Input* eps_xy_input  = element->GetInput(StrainRatexyEnum);     _assert_(eps_xy_input);
	Input* eps_yy_input  = element->GetInput(StrainRateyyEnum);     _assert_(eps_yy_input);
	Input*  n_input=element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
	Input* damage_input = NULL;
	Input* B_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
		B_input=element->GetInput(MaterialsRheologyBbarEnum); _assert_(B_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
		B_input=element->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	}

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		eps_xx_input->GetInputValue(&eps_xx,gauss);
		eps_xy_input->GetInputValue(&eps_xy,gauss);
		eps_yy_input->GetInputValue(&eps_yy,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		damage_input->GetInputValue(&damage,gauss);

		/*Calculate principal effective strain rates*/
		eps1=(eps_xx+eps_yy)/2.+sqrt(pow((eps_xx-eps_yy)/2.,2)+pow(eps_xy,2));
		eps2=(eps_xx+eps_yy)/2.-sqrt(pow((eps_xx-eps_yy)/2.,2)+pow(eps_xy,2));
		if(fabs(eps2)>fabs(eps1)){epstmp=eps2; eps2=eps1; eps1=epstmp;}

		/*Calculate effective strain rate and threshold strain rate*/
		epseff=1./sqrt(2.)*sqrt(eps1*eps1-eps1*eps2+eps2*eps2);
		eps0=pow(stress_threshold/B,n);

		if(epseff>eps0){
			f[i]=1.-pow(eps0/epseff,1./n)*exp(-(epseff-eps0)/(epsf-eps0))-damage;
		}
		else f[i]=0;

		/*Edits from MM*/
		if(f[i]>10.) f[i]=10.;
		if(f[i]<-10.) f[i]=-10.;
	}

	/*Add input*/
	element->AddInput(DamageFEnum,f,P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(f);
	delete gauss;
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInputPralong(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble c1,c2,c3,healing,stress_threshold;
	IssmDouble s_xx,s_xy,s_xz,s_yy,s_yz,s_zz,s1,s2,stmp;
	IssmDouble Chi,Psi,PosPsi,NegPsi;
	IssmDouble damage,tau_xx,tau_xy,tau_xz,tau_yy,tau_yz,tau_zz,stressMaxPrincipal;
	int equivstress,domaintype,dim;

	IssmDouble  k_sigma=0.0; /* k_sigma in exponential damage model*/
	int isexperiment=1;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*retrieve parameters:*/
	element->FindParam(&c1,DamageC1Enum);
	element->FindParam(&c2,DamageC2Enum);
	element->FindParam(&c3,DamageC3Enum);
	element->FindParam(&healing,DamageHealingEnum);
	element->FindParam(&stress_threshold,DamageStressThresholdEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Get problem dimension*/
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented");
	}
	/*Compute stress tensor and Stress Max Principal: */
	element->ComputeDeviatoricStressTensor();
	if(dim==3){
		/*Only works in 3d because the pressure is defined*/
		element->StressMaxPrincipalCreateInput();
	}
	/*retrieve what we need: */
	Input* tau_xx_input  = element->GetInput(DeviatoricStressxxEnum);     _assert_(tau_xx_input);
	Input* tau_xy_input  = element->GetInput(DeviatoricStressxyEnum);     _assert_(tau_xy_input);
	Input* tau_yy_input  = element->GetInput(DeviatoricStressyyEnum);     _assert_(tau_yy_input);
	Input* tau_xz_input  = NULL;
	Input* tau_yz_input  = NULL;
	Input* tau_zz_input  = NULL;
	Input* stressMaxPrincipal_input = NULL;
	if(dim==3){
		tau_xz_input  = element->GetInput(DeviatoricStressxzEnum);     _assert_(tau_xz_input);
		tau_yz_input  = element->GetInput(DeviatoricStressyzEnum);     _assert_(tau_yz_input);
		tau_zz_input  = element->GetInput(DeviatoricStresszzEnum);     _assert_(tau_zz_input);
		stressMaxPrincipal_input = element->GetInput(StressMaxPrincipalEnum); _assert_(stressMaxPrincipal_input);
	}
	Input* damage_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
	}

	/*retrieve the desired type of equivalent stress*/
	element->FindParam(&equivstress,DamageEquivStressEnum);

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		damage_input->GetInputValue(&damage,gauss);
		tau_xx_input->GetInputValue(&tau_xx,gauss);
		tau_xy_input->GetInputValue(&tau_xy,gauss);
		tau_yy_input->GetInputValue(&tau_yy,gauss);
		if(dim==3){
			tau_xz_input->GetInputValue(&tau_xz,gauss);
			tau_yz_input->GetInputValue(&tau_yz,gauss);
			tau_zz_input->GetInputValue(&tau_zz,gauss);
		}
		/*Calculate effective stress components*/
		s_xx=tau_xx/(1.-damage);
		s_xy=tau_xy/(1.-damage);
		s_yy=tau_yy/(1.-damage);
		if(dim==3){
			s_xz=tau_xz/(1.-damage);
			s_yz=tau_yz/(1.-damage);
			s_zz=tau_zz/(1.-damage);
		}

		/*Calculate k_simga value*/
		if(isexperiment==0){
			k_sigma = c3;
		}else if(isexperiment==1){
			/*FIXME: Hard coding. Parameters from Duddu et al. (2020)*/
			IssmDouble k1=-2.63; // unit: -
			IssmDouble k2=7.24*1e-6; // unit:  Pa-1
			IssmDouble tau_inv1=tau_xx + tau_yy;
			if (dim == 3){
				k_sigma += tau_zz;
			}

			k_sigma = k1 + k2*tau_inv1;
		}else{
			_error_("not implemented");
		}

		/*Calculate principal effective stresses*/
		if(dim==2){
			s1=(s_xx+s_yy)/2.+sqrt(pow((s_xx-s_yy)/2.,2)+pow(s_xy,2));
			s2=(s_xx+s_yy)/2.-sqrt(pow((s_xx-s_yy)/2.,2)+pow(s_xy,2));
			if(fabs(s2)>fabs(s1)){stmp=s2; s2=s1; s1=stmp;}

			if(equivstress==0){ /* von Mises */
				Chi=sqrt(s1*s1-s1*s2+s2*s2);
			}
			else if(equivstress==1){ /* max principal stress */
				Chi=s1;
			}
			else if(equivstress==2){ /* Hayhurst stress invariant */
				IssmDouble alpha=0.21;
				IssmDouble beta=0.63;
				Chi=alpha*s1 + beta*sqrt(s1*s1-s1*s2+s2*s2) + (1-alpha-beta)*(s1 + s2);
			}else{
				_error_("Not implemented");
			}
			Psi=Chi-stress_threshold;
			//NegPsi=max(-Chi,0.); /* healing only for compressive stresses */
			NegPsi=max(-Psi,0.); /* healing only for compressive stresses */
			PosPsi=max(Psi,0.);
			f[i]= c1*(pow(PosPsi,c2) - healing*pow(NegPsi,c2))*pow((1./(1.-damage)),k_sigma);
		}
		else{
			if(equivstress==1){/* max principal stress */
				stressMaxPrincipal_input->GetInputValue(&stressMaxPrincipal,gauss);
				Chi=stressMaxPrincipal/(1.-damage);
			}
			else if(equivstress==0){/* von Mises */
				Chi=sqrt(((s_xx-s_yy)*(s_xx-s_yy)+(s_yy-s_zz)*(s_yy-s_zz)+(s_zz-s_xx)*(s_zz-s_xx)+6.*(s_xy*s_xy+s_yz*s_yz+s_xz*s_xz))/2.);
			}
			else if(equivstress==2){ /* Hayhurst stress invariant */
				IssmDouble alpha=0.21;
				IssmDouble beta=0.63;
				Chi=alpha*s1 \
					 + beta*sqrt(((s_xx-s_yy)*(s_xx-s_yy)+(s_yy-s_zz)*(s_yy-s_zz)+(s_zz-s_xx)*(s_zz-s_xx)+6.*(s_xy*s_xy+s_yz*s_yz+s_xz*s_xz))/2.);
					 + (1-alpha-beta)*(s_xx + s_yy + s_zz);
			}
			Psi=Chi-stress_threshold;
			//NegPsi=max(-Chi,0.); /* healing only for compressive stresses */
			NegPsi=max(-Psi,0.); /* healing only for compressive stresses */
			PosPsi=max(Psi,0.);
			f[i]= c1*(pow(PosPsi,c2) - healing*pow(NegPsi,c2))*pow((1./(1.-damage)),k_sigma);
		}
	}
	/*Add input*/
	element->AddInput(DamageFEnum,f,P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(f);
	delete gauss;
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInputTest(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble c1,c2,c3,c4,healing,stress_threshold;
	IssmDouble s1, s2, s3;
	IssmDouble s_inv1;
	IssmDouble pressure;
	IssmDouble Chi,Psi,PosPsi,NegPsi;
	IssmDouble damage,sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,sigma_zz,stressMaxPrincipal;
	int equivstress,domaintype,dim;

	IssmDouble  k_sigma=0.0; /* k_sigma in exponential damage model*/
	int isdamage_exponent;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*retrieve parameters:*/
	element->FindParam(&c1,DamageC1Enum);
	element->FindParam(&c2,DamageC2Enum);
	element->FindParam(&c3,DamageC3Enum);
	element->FindParam(&c4,DamageC4Enum);
	element->FindParam(&healing,DamageHealingEnum);
	element->FindParam(&stress_threshold,DamageStressThresholdEnum);
	element->FindParam(&isdamage_exponent,DamageIsDamageExponentEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Get problem dimension*/
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain2DverticalEnum:   dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented");
	}

	/*retrieve what we need: */
	Input* damage_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
	}

	Input* stress_equivalent_input = element->GetInput(DamageEffectiveStressEquivalentEnum); _assert_(stress_equivalent_input);
	/* stress */
	Input* sigma1_input = element->GetInput(StressTensorPrincipalValue1Enum); _assert_(sigma1_input);
	Input* sigma2_input = element->GetInput(StressTensorPrincipalValue2Enum); _assert_(sigma2_input);
	Input* sigma3_input = element->GetInput(StressTensorPrincipalValue3Enum); _assert_(sigma3_input);
	/* Effective stress */
	Input* effective_sigma1_input = element->GetInput(DamageEffectiveStressPrincipalValue1Enum); _assert_(effective_sigma1_input);
	Input* effective_sigma2_input = element->GetInput(DamageEffectiveStressPrincipalValue2Enum); _assert_(effective_sigma2_input);
	Input* effective_sigma3_input = element->GetInput(DamageEffectiveStressPrincipalValue3Enum); _assert_(effective_sigma3_input);

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		/*Get all variables*/
		damage_input->GetInputValue(&damage,gauss);
		stress_equivalent_input->GetInputValue(&Chi,gauss);

		/*Calculate k_simga value*/
		if(isdamage_exponent==0){
			k_sigma = c3;
		}else if(isdamage_exponent==1){
			/* Eq 23. Pralong et al. (2005) */
			IssmDouble k1=3.75e-3; // [Pa-1] Table 3 in Pralong et al. (2005)

			effective_sigma1_input->GetInputValue(&s1,gauss);
			effective_sigma2_input->GetInputValue(&s2,gauss);
			effective_sigma3_input->GetInputValue(&s3,gauss);
			s_inv1 = s1 + s2 + s3;

			k_sigma = k1 * ( sqrt(max(s_inv1,0.0)) - healing * sqrt(max(-s_inv1,0.0)) );
		}else if(isdamage_exponent==2){
			/*FIXME: Hard coding. Parameters from Duddu et al. (2013): doi:10.1007/s00466-012-0778-7 */
			_error_("Error: not supporte expornent in damage (k_sigma)");
			IssmDouble k1=-2.63;
			IssmDouble k2=7.24*1e-6; // unit: Pa-1
			IssmDouble k3=12.36;
			IssmDouble k4=-3.83; // unit: Pa-1 

			sigma1_input->GetInputValue(&s1,gauss);
			sigma2_input->GetInputValue(&s2,gauss);
			sigma3_input->GetInputValue(&s3,gauss);
			s_inv1 = s1+s2+s3;

			/* Eq 18 in Duddu 2013 */
			if (s_inv1 >= 0) k_sigma = k1 + k2*abs(s_inv1); // tension
			else if(s_inv1 <0) k_sigma = -(k3 + k3*abs(s_inv1)); // compression
			else _error_("Error: out of bounds");
		}else if(isdamage_exponent==3){
			/*FIXME: Hard coding. Parameters from Duddu et al. (2020)*/
			IssmDouble k1=-2.63; // unit: -
			IssmDouble k2=7.24*1e-6; // unit: Pa-1

			sigma1_input->GetInputValue(&s1,gauss);
			sigma2_input->GetInputValue(&s2,gauss);
			sigma3_input->GetInputValue(&s3,gauss);
			s_inv1 = s1+s2+s3;

			k_sigma = k1 + k2*s_inv1;
		}else{
			_error_("not implemented");
		}

		/*Calculate principal effective stresses*/
		if(dim==2){
			Psi=Chi-stress_threshold;
			PosPsi=max(Psi,0.);
			NegPsi=max(-Psi,0.); /* healing only for compressive stresses */
			f[i]= c1*(pow(PosPsi,c2) - healing*pow(NegPsi,c2))*pow((1./(1.-damage)),k_sigma);
		}
		else{
			Psi=Chi-stress_threshold;
			PosPsi=max(Psi,0.);
			NegPsi=max(-Psi,0.); /* healing only for compressive stresses */
			f[i]= c1*(pow(PosPsi,c2) - healing*pow(NegPsi,c2))*pow((1./(1.-damage)),k_sigma);
		}
	}
	/*Add input*/
	element->AddInput(DamageFEnum,f,P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(f);
	delete gauss;
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInputLinear(Element* element){/*{{{*/
	/*
	 Damage evolution as follows:

	 f = B * max(Chi,0);
	 */

	/*Intermediaries */
	IssmDouble c1,healing,stress_threshold;
	IssmDouble s_xx,s_xy,s_xz,s_yy,s_yz,s_zz,s1,s2,s3,stmp;
	IssmDouble J2s,Chi,Psi,PosPsi,NegPsi;
	IssmDouble damage,sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,sigma_zz,stressMaxPrincipal;
	int equivstress,domaintype,dim;

	IssmDouble  k_sigma=0.0; /* k_sigma in exponential damage model*/
	int isexperiment=1;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*retrieve parameters:*/
	element->FindParam(&c1,DamageC1Enum);
	element->FindParam(&healing,DamageHealingEnum);
	element->FindParam(&stress_threshold,DamageStressThresholdEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Get problem dimension*/
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented");
	}
	/*Compute stress tensor and Stress Max Principal: */
	element->ComputeStressTensor();
	if(dim==3){
		/*Only works in 3d because the pressure is defined*/
		element->StressMaxPrincipalCreateInput();
	}
	/*retrieve what we need: */
	Input* sigma_xx_input  = element->GetInput(StressTensorxxEnum);     _assert_(sigma_xx_input);
	Input* sigma_xy_input  = element->GetInput(StressTensorxyEnum);     _assert_(sigma_xy_input);
	Input* sigma_yy_input  = element->GetInput(StressTensoryyEnum);     _assert_(sigma_yy_input);
	Input* sigma_xz_input  = NULL;
	Input* sigma_yz_input  = NULL;
	Input* sigma_zz_input  = NULL;
	Input* stressMaxPrincipal_input = NULL;
	if(dim==3){
		sigma_xz_input  = element->GetInput(StressTensorxzEnum);     _assert_(sigma_xz_input);
		sigma_yz_input  = element->GetInput(StressTensoryzEnum);     _assert_(sigma_yz_input);
		sigma_zz_input  = element->GetInput(StressTensorzzEnum);     _assert_(sigma_zz_input);
		stressMaxPrincipal_input = element->GetInput(StressMaxPrincipalEnum); _assert_(stressMaxPrincipal_input);
	}
	Input* damage_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
	}

	/*retrieve the desired type of equivalent stress*/
	element->FindParam(&equivstress,DamageEquivStressEnum);

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		damage_input->GetInputValue(&damage,gauss);
		sigma_xx_input->GetInputValue(&sigma_xx,gauss);
		sigma_xy_input->GetInputValue(&sigma_xy,gauss);
		sigma_yy_input->GetInputValue(&sigma_yy,gauss);
		if(dim==3){
			sigma_xz_input->GetInputValue(&sigma_xz,gauss);
			sigma_yz_input->GetInputValue(&sigma_yz,gauss);
			sigma_zz_input->GetInputValue(&sigma_zz,gauss);
		}
		/*Calculate effective stress components*/
		s_xx=sigma_xx/(1.-damage);
		s_xy=sigma_xy/(1.-damage);
		s_yy=sigma_yy/(1.-damage);
		if(dim==3){
			s_xz=sigma_xz/(1.-damage);
			s_yz=sigma_yz/(1.-damage);
			s_zz=sigma_zz/(1.-damage);
		}

		/*Calculate principal effective stresses*/
		if(dim==2){
			/* Compute principal effective stresses */
			Matrix2x2Eigen(&s1,&s2,NULL,NULL,s_xx,s_xy,s_yy);

			if(s2>s1){stmp=s2; s2=s1; s1=stmp;}

			if(equivstress==0){ /* von Mises */
				Chi=sqrt(s1*s1-s1*s2+s2*s2);
			}
			else if(equivstress==1){ /* max principal stress */
				Chi=s1;
			}
			else if(equivstress==2){ /* Hayhurst stress invariant */
				IssmDouble alpha=0.21;
				IssmDouble beta=0.63;
				Chi=alpha*s1 + beta*sqrt(s1*s1-s1*s2+s2*s2) + (1-alpha-beta)*(s1 + s2);
			}else{
				_error_("Not implemented");
			}
			Psi=Chi-stress_threshold;
			PosPsi=max(Psi,0.);
			f[i]= c1*PosPsi;
		}
		else{
			_error_("Not implemented yet.");
			if(equivstress==1){/* max principal stress */
				stressMaxPrincipal_input->GetInputValue(&stressMaxPrincipal,gauss);
				Chi=stressMaxPrincipal/(1.-damage);
			}
			else if(equivstress==0){/* von Mises */
				Chi=sqrt(((s_xx-s_yy)*(s_xx-s_yy)+(s_yy-s_zz)*(s_yy-s_zz)+(s_zz-s_xx)*(s_zz-s_xx)+6.*(s_xy*s_xy+s_yz*s_yz+s_xz*s_xz))/2.);
			}
			else if(equivstress==2){ /* Hayhurst stress invariant */
				IssmDouble alpha=0.21;
				IssmDouble beta=0.63;
				Chi=alpha*s1 \
					 + beta*sqrt(((s_xx-s_yy)*(s_xx-s_yy)+(s_yy-s_zz)*(s_yy-s_zz)+(s_zz-s_xx)*(s_zz-s_xx)+6.*(s_xy*s_xy+s_yz*s_yz+s_xz*s_xz))/2.);
					 + (1-alpha-beta)*(s_xx + s_yy + s_zz);
			}
			Psi=Chi-stress_threshold;
			NegPsi=max(-Psi,0.); /* healing only for compressive stresses */
			PosPsi=max(Psi,0.);
			f[i] = c1*PosPsi;
		}
	}
	/*Add input*/
	element->AddInput(DamageFEnum,f,P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(f);
	delete gauss;
}/*}}}*/
ElementVector* DamageEvolutionAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* DamageEvolutionAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* DamageEvolutionAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;
	/*Intermediaries*/
	int         domaintype,dim;
	int         stabilization;
	IssmDouble  Jdet,dt,D_scalar,h,hx,hy,hz;
	IssmDouble  vel,vx,vy,vz,dvxdx,dvydy,dvzdz; //dvx[3],dvy[3],dvz[3];
	IssmDouble  u,v,w,um,vm,wm;
	IssmDouble  du[3],dv[3],dw[3],dum[3],dvm[3],dwm[3];
	IssmDouble *xyz_list  = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("Not implemented yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,DamageStabilizationEnum);
	Input* vx_input = element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = NULL;
	Input* vxm_input= element->GetInput(VxMeshEnum); _assert_(vxm_input);
	Input* vym_input= element->GetInput(VyMeshEnum); _assert_(vym_input);
	Input* vzm_input= NULL;
	if(dim==3){
		vz_input =element->GetInput(VzEnum);     _assert_(vz_input);
		vzm_input=element->GetInput(VzMeshEnum); _assert_(vzm_input);
	}

	if(dim==2) h=element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		vx_input->GetInputValue(&u,gauss); vxm_input->GetInputValue(&um,gauss); vx=u-um;
		vy_input->GetInputValue(&v,gauss); vym_input->GetInputValue(&vm,gauss); vy=v-vm;
		vx_input->GetInputDerivativeValue(&du[0],xyz_list,gauss); vxm_input->GetInputDerivativeValue(&dum[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dv[0],xyz_list,gauss); vym_input->GetInputDerivativeValue(&dvm[0],xyz_list,gauss);
		if(dim==3){
			vz_input->GetInputValue(&w,gauss); vz_input->GetInputValue(&wm,gauss); vz=w-wm;
			vz_input->GetInputDerivativeValue(&dw[0],xyz_list,gauss); vzm_input->GetInputDerivativeValue(&dwm[0],xyz_list,gauss);
		}

		/*Transient term*/
		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j];

		/*Advection: */
		dvxdx=du[0]-dum[0];
		dvydy=dv[1]-dvm[1];
		D_scalar=dt*gauss->weight*Jdet;
		if(dim==2){
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*\phi_i \phi_j \nabla\cdot v*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy);
					/*\phi_i v\cdot\nabla\phi_j*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
				}
			}
		}
		else{/*3D*/
			_assert_(dim==3);
			dvzdz=dw[2]-dwm[2];
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*\phi_i \phi_j \nabla\cdot v*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy+dvzdz);
					/*\phi_i v\cdot\nabla\phi_j*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j] + vz*dbasis[2*numnodes+j]);
				}
			}
		}

		/*Streamline upwinding: */
		if(stabilization==2){
			if(dim==3){
				vel=sqrt(vx*vx+vy*vy+vz*vz)+1.e-8;
				D[0*dim+0]=h/(2.0*vel)*vx*vx;
				D[1*dim+0]=h/(2.0*vel)*vy*vx;
				D[2*dim+0]=h/(2.0*vel)*vz*vx;
				D[0*dim+1]=h/(2.0*vel)*vx*vy;
				D[1*dim+1]=h/(2.0*vel)*vy*vy;
				D[2*dim+1]=h/(2.0*vel)*vy*vz;
				D[0*dim+2]=h/(2.0*vel)*vx*vz;
				D[1*dim+2]=h/(2.0*vel)*vy*vz;
				D[2*dim+2]=h/(2.0*vel)*vz*vz;
			}
			else if(dim=2){
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				D[0*dim+0]=h/(2.0*vel)*vx*vx;
				D[1*dim+0]=h/(2.0*vel)*vy*vx;
				D[0*dim+1]=h/(2.0*vel)*vx*vy;
				D[1*dim+1]=h/(2.0*vel)*vy*vy;
			}
		}
		/*Artificial diffusivity: */
		else if(stabilization==1){
			if(dim==2){
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);
				D[0*dim+0]=h/2.0*fabs(vx);
				D[1*dim+1]=h/2.0*fabs(vy);
			}
			else if(dim==3){
				element->ElementSizes(&hx,&hy,&hz);
				vel=sqrt(vx*vx + vy*vy + vz*vz)+1.e-14;
				h=sqrt( pow(hx*vx/vel,2) + pow(hy*vy/vel,2) + pow(hz*vz/vel,2));
				D[0*dim+0]=h/(2.*vel)*fabs(vx*vx);  D[0*dim+1]=h/(2.*vel)*fabs(vx*vy); D[0*dim+2]=h/(2.*vel)*fabs(vx*vz);
				D[1*dim+0]=h/(2.*vel)*fabs(vy*vx);  D[1*dim+1]=h/(2.*vel)*fabs(vy*vy); D[1*dim+2]=h/(2.*vel)*fabs(vy*vz);
				D[2*dim+0]=h/(2.*vel)*fabs(vz*vx);  D[2*dim+1]=h/(2.*vel)*fabs(vz*vy); D[2*dim+2]=h/(2.*vel)*fabs(vz*vz);
			}
		}
		if(stabilization==1 || stabilization==2){
			if(dim==2){
				D[0*dim+0]=D_scalar*D[0*dim+0];
				D[1*dim+0]=D_scalar*D[1*dim+0];
				D[0*dim+1]=D_scalar*D[0*dim+1];
				D[1*dim+1]=D_scalar*D[1*dim+1];
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += (
									dbasis[0*numnodes+i] *(D[0*dim+0]*dbasis[0*numnodes+j] + D[0*dim+1]*dbasis[1*numnodes+j]) +
									dbasis[1*numnodes+i] *(D[1*dim+0]*dbasis[0*numnodes+j] + D[1*dim+1]*dbasis[1*numnodes+j])
									);
					}
				}
			}
			else if(dim==3){
				D[0*dim+0]=D_scalar*D[0*dim+0];
				D[1*dim+0]=D_scalar*D[1*dim+0];
				D[2*dim+0]=D_scalar*D[2*dim+0];
				D[0*dim+1]=D_scalar*D[0*dim+1];
				D[1*dim+1]=D_scalar*D[1*dim+1];
				D[2*dim+1]=D_scalar*D[2*dim+1];
				D[0*dim+2]=D_scalar*D[0*dim+2];
				D[1*dim+2]=D_scalar*D[1*dim+2];
				D[2*dim+2]=D_scalar*D[2*dim+2];
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += (
									dbasis[0*numnodes+i] *(D[0*dim+0]*dbasis[0*numnodes+j] + D[0*dim+1]*dbasis[1*numnodes+j] + D[0*dim+2]*dbasis[2*numnodes+j]) +
									dbasis[1*numnodes+i] *(D[1*dim+0]*dbasis[0*numnodes+j] + D[1*dim+1]*dbasis[1*numnodes+j] + D[1*dim+2]*dbasis[2*numnodes+j]) +
									dbasis[2*numnodes+i] *(D[2*dim+0]*dbasis[0*numnodes+j] + D[2*dim+1]*dbasis[1*numnodes+j] + D[2*dim+2]*dbasis[2*numnodes+j])
									);
					}
				}
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(D);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* DamageEvolutionAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype,damagelaw;
	IssmDouble  Jdet,dt;
	IssmDouble  f,damage;
	IssmDouble* xyz_list = NULL;
	/*Get element*/
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&damagelaw,DamageLawEnum);
	switch(damagelaw){
		case 0:
			this->CreateDamageFInput(element);
			break;
		case 1:
			this->CreateDamageFInputPralong(element);
			break;
		case 2:
			this->CreateDamageFInputExp(element);
			break;
		case 3:
			this->CreateDamageFInputArctan(element);
			break;
		case 4:
			this->CreateDamageFInputTest(element);
			break;
		case 5:
			this->CreateDamageFInputLinear(element);
			break;
		default:
			_error_("not implemented yet");
	}

	Input* damaged_input = NULL;
	Input* damagef_input = element->GetInput(DamageFEnum); _assert_(damagef_input);
	if(domaintype==Domain2DhorizontalEnum){
		damaged_input = element->GetInput(DamageDbarEnum); _assert_(damaged_input);
	}
	else{
		damaged_input = element->GetInput(DamageDEnum); _assert_(damaged_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		damaged_input->GetInputValue(&damage,gauss);
		damagef_input->GetInputValue(&f,gauss);

		IssmDouble factor = Jdet*gauss->weight*(damage+dt*f);
		for(int i=0;i<numnodes;i++){
			pe->values[i]+=factor*basis[i];
		}
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           DamageEvolutionAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,DamageDbarEnum);
}/*}}}*/
void           DamageEvolutionAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           DamageEvolutionAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	int *doflist = NULL;
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* newdamage = xNew<IssmDouble>(numnodes);

	/*Get user-supplied max_damage: */
	IssmDouble max_damage = element->FindParam(DamageMaxDamageEnum);
	IssmDouble critical_damage = element->FindParam(DamageCriticalDamageEnum);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		newdamage[i]=solution[doflist[i]];
		/*Check solution*/
		if(xIsNan<IssmDouble>(newdamage[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(newdamage[i])) _error_("Inf found in solution vector");
		/*Enforce D < max_damage and D > 0 */
		if(newdamage[i]>max_damage) newdamage[i]=max_damage;
		else if(newdamage[i]<0.)    newdamage[i]=0.;

		/*Enforce max_damage if D > critical_damage is satisfied (Huth 2021, doi:10.1029/2020MS002292) */
		if(newdamage[i]>critical_damage) newdamage[i]=max_damage;
	}

	/*Get all inputs and parameters*/
	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype==Domain2DhorizontalEnum){
		element->AddInput(DamageDbarEnum,newdamage,element->GetElementType());
	}
	else{
		element->AddInput(DamageDEnum,newdamage,element->GetElementType());
	}

	/*Free resources:*/
	xDelete<IssmDouble>(newdamage);
	xDelete<int>(doflist);
}/*}}}*/
void           DamageEvolutionAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/

/*Flux Correction Transport*/
ElementMatrix* DamageEvolutionAnalysis::CreateFctKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble Jdet;
	IssmDouble vx,vy;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int dim      = 2;

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vxaverage_input=element->GetInput(VxEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=element->GetInput(VyEnum); _assert_(vyaverage_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);

		IssmDouble factor = -gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/*\phi_i v\cdot\nabla\phi_j*/
				Ke->values[i*numnodes+j] += factor*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
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
ElementMatrix* DamageEvolutionAnalysis::CreateMassMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Me     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		D=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Me->values[i*numnodes+j] += D*basis[i]*basis[j];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Me;
}/*}}}*/
ElementVector* DamageEvolutionAnalysis::CreateFctPVector(Element* element){/*{{{*/

	return this->CreatePVector(element);

}/*}}}*/
void           DamageEvolutionAnalysis::FctKMatrix(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,FemModel* femmodel){/*{{{*/

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
void           DamageEvolutionAnalysis::FctPVector(Vector<IssmDouble>** ppf,FemModel* femmodel){/*{{{*/

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
void           DamageEvolutionAnalysis::MassMatrix(Matrix<IssmDouble>** pMff,FemModel* femmodel){/*{{{*/

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

/*Others*/
void DamageEvolutionAnalysis::ComputeEffectiveStress(FemModel* femmodel){/*{{{*/

	/* Compute effective stress equivalent values at once */
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		ComputeEffectiveStress(element);
	}

}/*}}}*/
void DamageEvolutionAnalysis::ComputeEffectiveStress(Element* element){/*{{{*/
	/* Compute the equivalent stress for determining the damage evolution and store values relavant to effective stress tensor: principal and invariants
	*/
	
	/* Intermediaries */
	int domaintype, dim;
	element->FindParam(&domaintype,DomainTypeEnum);

	/* Get problem dimension */
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain2DverticalEnum:   dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented yet");
	}
	
	/* Intermediaries */
	IssmDouble tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz; /* Deviatoric stress tensor components [Pa] */
	IssmDouble sigma_xx, sigma_xy, sigma_xz, sigma_yy, sigma_yz, sigma_zz; /* Cauchy stress tensor components [Pa] */
	IssmDouble s1, s2, s3; /* principal effective Cauchy stress [Pa] */
	IssmDouble t1, t2, t3; /* principal effective devitaoric stress [Pa] */
	IssmDouble stmp;
	IssmDouble D; /* Damage value [unit: 1]*/
	IssmDouble z; /* elevation [m] */

	IssmDouble rho_w, g;
	IssmDouble alpha, beta; /* parameter for Hayhurst */
	IssmDouble mu; /* parameter for Coulomb */

	/*    Principal stress for each vertex */
	IssmDouble *sigma1 = NULL;
	IssmDouble *sigma2 = NULL;
	IssmDouble *sigma3 = NULL;
	IssmDouble *effective_sigma1 = NULL;
	IssmDouble *effective_sigma2 = NULL;
	IssmDouble *effective_sigma3 = NULL;
	IssmDouble *effective_tau1 = NULL;
	IssmDouble *effective_tau2 = NULL;
	IssmDouble *effective_tau3 = NULL;
	IssmDouble *xyz_list = NULL;
	
	/*    ice flow equation */
	bool isSSA, isL1L2, isMOLHO, isHO, isFS;
	
	/* Pressure term
	P  : pressure for stress tensor [Pa]
	Pi : Overburden pressure [Pa] 
	Pw : Water pressure [Pa]
	*/
	IssmDouble P, Pi, Pw, sealevel;
	bool isPeff; /* if true, we compute the effective stress, else we compute the Cauchy stress */
	int iseffective_stress; /* effective stress. 0: pralong2005, 1: keller 2014 */

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();

	sigma1  = xNewZeroInit<IssmDouble>(numnodes);
	sigma2  = xNewZeroInit<IssmDouble>(numnodes);
	sigma3  = xNewZeroInit<IssmDouble>(numnodes);

	effective_sigma1  = xNewZeroInit<IssmDouble>(numnodes);
	effective_sigma2  = xNewZeroInit<IssmDouble>(numnodes);
	effective_sigma3  = xNewZeroInit<IssmDouble>(numnodes);

	effective_tau1  = xNewZeroInit<IssmDouble>(numnodes);
	effective_tau2  = xNewZeroInit<IssmDouble>(numnodes);
	effective_tau3  = xNewZeroInit<IssmDouble>(numnodes);

	/* Retrieve constants */
	element->FindParam(&rho_w,MaterialsRhoSeawaterEnum);
	element->FindParam(&g,ConstantsGEnum);
	element->FindParam(&isPeff,DamageIsPeffEnum);
	element->FindParam(&iseffective_stress,DamageIsEffectiveStressEnum);

	element->FindParam(&isSSA,FlowequationIsSSAEnum);
	element->FindParam(&isL1L2,FlowequationIsL1L2Enum);
	element->FindParam(&isMOLHO,FlowequationIsMOLHOEnum);
	element->FindParam(&isHO,FlowequationIsHOEnum);
	element->FindParam(&isFS,FlowequationIsFSEnum);

	/* Retrieve all input and parameters */
	element->GetVerticesCoordinates(&xyz_list);

	/* Precompute deviatoric stress tensor*/
	element->ComputeDeviatoricStressTensor();

	/* Retrieve what we need: */
	Input* tau_xx_input  = element->GetInput(DeviatoricStressxxEnum); _assert_(tau_xx_input);
	Input* tau_xy_input  = element->GetInput(DeviatoricStressxyEnum); _assert_(tau_xy_input);
	Input* tau_yy_input  = element->GetInput(DeviatoricStressyyEnum); _assert_(tau_yy_input);
	Input* tau_xz_input  = NULL;
	Input* tau_yz_input  = NULL;
	Input* tau_zz_input  = NULL;
	if(dim==3){
		tau_xz_input  = element->GetInput(DeviatoricStressxzEnum);     _assert_(tau_xz_input);
		tau_yz_input  = element->GetInput(DeviatoricStressyzEnum);     _assert_(tau_yz_input);
		tau_zz_input  = element->GetInput(DeviatoricStresszzEnum);     _assert_(tau_zz_input);
	}
	/* NOTE: pressure in model would be overburden pressure for SSA and HO */
	Input* pressure_input = element->GetInput(PressureEnum); _assert_(pressure_input);
	Input* sealevel_input = element->GetInput(SealevelEnum); _assert_(sealevel_input);

	Input* damage_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){ damage_input = element->GetInput(DamageDbarEnum); _assert_(damage_input); }
	else if(domaintype==Domain2DverticalEnum){ damage_input = element->GetInput(DamageDEnum); _assert_(damage_input); }
	else{ damage_input = element->GetInput(DamageDEnum); _assert_(damage_input); }

	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		/* Get variables */
		damage_input->GetInputValue(&D,gauss);
		sealevel_input->GetInputValue(&sealevel,gauss);

		/* Get deviatoric stress */
		tau_xx_input->GetInputValue(&tau_xx,gauss);
		tau_xy_input->GetInputValue(&tau_xy,gauss);
		tau_yy_input->GetInputValue(&tau_yy,gauss);
		if(dim==3){
			tau_xz_input->GetInputValue(&tau_xz,gauss);
			tau_yz_input->GetInputValue(&tau_yz,gauss);
			tau_zz_input->GetInputValue(&tau_zz,gauss);
		}
		
		/* Compute pressure in Cauchy stress tensor*/
		if (isSSA | isHO){ /* {{{ */
			pressure_input->GetInputValue(&Pi,gauss);
			if(dim==2){
				_error_("Error: not supported.");
				P = 0.0; /* pressure at surface*/
			}else if(dim==3){
				z = xyz_list[i*3+2];
				Pw = 0.0;
				if (isPeff){
					if (z > sealevel){ /* explicitly set Pw value */
						Pw = 0.0;
					}else if(z <= sealevel){
						Pw = rho_w*g*(sealevel-z);
					}
				}
				/* Eq. 15 in Huth et al. (2021) / SSA and HO assume the hydrostatic assumption (Greve and Blatter, 2009) */
				if (iseffective_stress == 0) P = Pi - (tau_xx + tau_yy) - Pw;
				else if (iseffective_stress == 1) P = Pi - (tau_xx + tau_yy)/(1-D) - Pw;
			}
		}/* }}} */
		else if(isFS){ /* {{{ */
			pressure_input->GetInputValue(&P,gauss);
		}/* }}} */

		/* Compute effective Cauchy stress tensor based on deviatoric stress */
		if (iseffective_stress == 0){ /* Pralong 2005 */
			sigma_xx = (tau_xx-P)/(1-D);
			sigma_xy = (tau_xy)/(1-D);
			sigma_yy = (tau_yy-P)/(1-D);
			if(dim==3){
				sigma_xz = tau_xz/(1-D);
				sigma_yz = tau_yz/(1-D);
				//if(isSSA | isHO) sigma_zz = (tau_zz-P)/(1-D);
				if(isSSA | isHO) sigma_zz = (-tau_xx - tau_yy - P)/(1-D);
				else if(isFS) sigma_zz = (tau_zz-P)/(1-D);
				else _error_("Error: Not supported.");
			}
		}else if (iseffective_stress == 1){ /* Keller 2014 */
			sigma_xx = tau_xx/(1-D) - P;
			sigma_xy = tau_xy/(1-D);
			sigma_yy = tau_yy/(1-D) - P;
			if(dim==3){
				sigma_xz = tau_xz/(1-D);
				sigma_yz = tau_yz/(1-D);
				//if(isSSA | isHO) sigma_zz = tau_zz/(1-D) - P;
				if(isSSA | isHO) sigma_zz = (-tau_xx - tau_yy)/(1-D) - P;
				else if(isFS) sigma_zz = tau_zz/(1-D) - P;
				else _error_("Error: Not supported.");
			}
		}

		/* Compute principal stress values without damage */
		if(dim==2){/* {{{ */
			_error_("not supported yet");
		}/* }}} */
		else if(dim==3){/* {{{ */
			if(isSSA){ /* Hydrostatic assumption : neglect sigma_zx, sigma_zy */
				/* directly compute stress tensor */
				Matrix2x2Eigen(&s1,&s2,NULL,NULL,2*tau_xx+tau_yy-Pi+Pw, tau_xy, tau_xx+2*tau_yy-Pi+Pw);
				s3 = -Pi+Pw;
			}else if(isHO){
				Matrix3x3Eigen(&s1,&s2,&s3,NULL,NULL,NULL,
					2*tau_xx+tau_yy-Pi+Pw, tau_xy               , tau_xz,
					tau_xy               , tau_xx+2*tau_yy-Pi+Pw, tau_yz,
					0.0                  , 0.0                  , -Pi+Pw);
			}else if(isFS){
				Matrix3x3Eigen(&s1,&s2,&s3,NULL,NULL,NULL,
					tau_xx-P, tau_xy  , tau_xz,
					tau_xy  , tau_yy-P, tau_yz,
					tau_xz  , tau_yz  , tau_zz-P);
			}else{
				_error_("not implemented yet");
			}

			/* Ordering large value (descending order) ... s1 > s2 > s3 */
			if(s1<s2) swap(s1,s2);
			if(s1<s3) swap(s1,s3);
			if(s2<s3) swap(s2,s3);
			_assert_( (s1>=s2) & (s2 >= s3) );
		}/* }}} */

		/* Now, assign principal stresses */
		sigma1[i] = s1;
		sigma2[i] = s2;
		sigma3[i] = s3;
		
		/* Compute effective deviatoric stress */
		tau_xx = tau_xx/(1-D);
		tau_xy = tau_xy/(1-D);
		tau_yy = tau_yy/(1-D);
		if(dim==3){
			tau_xz = tau_xz/(1-D);
			tau_yz = tau_yz/(1-D);
			tau_zz = tau_zz/(1-D);
		}

		/* Compute effective principal stress values */
		if(dim==2){ /* {{{ */
			_error_("Error: not supported?!");
			/* Compute principal effective stresses */
			Matrix2x2Eigen(&s1,&s2,NULL,NULL,sigma_xx,sigma_xy,sigma_yy);
			s3 = sigma_zz;
			Matrix2x2Eigen(&t1,&t2,NULL,NULL,tau_xx,tau_xy,tau_yy);
			t3 = tau_zz;
			/* }}} */
		}else if(dim==3){ /* {{{ */
			/* Compute principal effective stresses*/
			/* FIXME: Now, stress equivalent would be only computed for SSA and HO, not FS.*/
			if(isSSA){ /* Hydrostatic assumption : neglect sigma_zx, sigma_zy */
				Matrix2x2Eigen(&s1,&s2,NULL,NULL,sigma_xx,sigma_xy,sigma_yy);
				s3 = sigma_zz;
				//Matrix3x3Eigen(&s1,&s2,&s3,NULL,NULL,NULL,
				//	sigma_xx,sigma_xy,sigma_xz,
				//	sigma_xy,sigma_yy,sigma_yz,
				//	0.0, 0.0, sigma_zz);
				Matrix2x2Eigen(&t1,&t2,NULL,NULL,tau_xx,tau_xy,tau_yy);
				t3 = tau_zz;
				//Matrix3x3Eigen(&t1,&t2,&t3,NULL,NULL,NULL,
				//	tau_xx,tau_xy,0.0,
				//	tau_xy,tau_yy,0.0,
				//	0.0, 0.0, tau_zz);
			}else if(isHO){
				Matrix3x3Eigen(&s1,&s2,&s3,NULL,NULL,NULL,
					sigma_xx,sigma_xy,sigma_xz,
					sigma_xy,sigma_yy,sigma_yz,
					0.0, 0.0, sigma_zz);
				Matrix3x3Eigen(&t1,&t2,&t3,NULL,NULL,NULL,
					tau_xx,tau_xy,tau_xz,
					tau_xy,tau_yy,tau_yz,
					0.0, 0.0, tau_zz);
			}else if(isFS){
				_error_("Not implemented for FS equation.");
			}else{
				_error_("not implemented yet");
			}

			/* Ordering large value (descending order) ... s1 > s2 > s3 */
			if(s1<s2) swap(s1,s2);
			if(s1<s3) swap(s1,s3);
			if(s2<s3) swap(s2,s3);
			_assert_( (s1>=s2) & (s2 >= s3) );

			if(t1<t2) swap(t1,t2);
			if(t1<t3) swap(t1,t3);
			if(t2<t3) swap(t2,t3);
			_assert_( (t1>=t2) & (t2 >= t3) );

		} /* }}} */
		
		/* Now, assign effective principal stresses */
		effective_sigma1[i] = s1;
		effective_sigma2[i] = s2;
		effective_sigma3[i] = s3;

		effective_tau1[i] = t1;
		effective_tau2[i] = t2;
		effective_tau3[i] = t3;
	}

	/* Assign values */
	element->AddInput(StressTensorPrincipalValue1Enum,sigma1,P1DGEnum);
	element->AddInput(StressTensorPrincipalValue2Enum,sigma2,P1DGEnum);
	element->AddInput(StressTensorPrincipalValue3Enum,sigma3,P1DGEnum);

	element->AddInput(DamageEffectiveStressPrincipalValue1Enum,effective_sigma1,P1DGEnum);
	element->AddInput(DamageEffectiveStressPrincipalValue2Enum,effective_sigma2,P1DGEnum);
	element->AddInput(DamageEffectiveStressPrincipalValue3Enum,effective_sigma3,P1DGEnum);

	element->AddInput(DamageEffectiveDeviatoricStressPrincipalValue1Enum,effective_tau1,P1DGEnum);
	element->AddInput(DamageEffectiveDeviatoricStressPrincipalValue2Enum,effective_tau2,P1DGEnum);
	element->AddInput(DamageEffectiveDeviatoricStressPrincipalValue3Enum,effective_tau3,P1DGEnum);

	/* Clear memory */
	xDelete<IssmDouble>(sigma1);
	xDelete<IssmDouble>(sigma2);
	xDelete<IssmDouble>(sigma3);

	xDelete<IssmDouble>(effective_sigma1);
	xDelete<IssmDouble>(effective_sigma2);
	xDelete<IssmDouble>(effective_sigma3);
	xDelete<IssmDouble>(effective_tau1);
	xDelete<IssmDouble>(effective_tau2);
	xDelete<IssmDouble>(effective_tau3);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
}/*}}}*/
void DamageEvolutionAnalysis::ComputeStressEquivalent(FemModel* femmodel){/*{{{*/
	
	/* Compute all stress equivalent values at once */
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		ComputeStressEquivalent(element);
	}

}/*}}}*/
void DamageEvolutionAnalysis::ComputeStressEquivalent(Element* element){/*{{{*/
	/* Compute the equivalent stress for determining the damage evolution and store values relavant to effective stress tensor: principal and invariants
	*/
	
	/* Intermediaries */
	int domaintype, dim;
	element->FindParam(&domaintype,DomainTypeEnum);

	/* Get problem dimension */
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain2DverticalEnum:   dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented yet");
	}
	
	/* Intermediaries */
	IssmDouble s1, s2, s3; /* principal effective Cauchy stress [Pa] */
	IssmDouble t1, t2, t3; /* principal effective devitaoric stress [Pa] */
	IssmDouble *sigma_equiv;
	IssmDouble P; /* pressure */
	int isequivstress; /* stress equivalent stress */
	int iseffective_stress; /* scheme for effective stress */
	double alpha, beta, mu;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	sigma_equiv = xNewZeroInit<IssmDouble>(numnodes);

	/* Retrieve constants */
	element->FindParam(&isequivstress,DamageEquivStressEnum);
	element->FindParam(&iseffective_stress,DamageIsEffectiveStressEnum);
	element->FindParam(&alpha,DamageEquivStressAlphaEnum);
	element->FindParam(&beta,DamageEquivStressBetaEnum);
	element->FindParam(&mu,DamageEquivStressMuEnum);

	/* Retrieve what we need: */
	Input* sigma1_input = element->GetInput(DamageEffectiveStressPrincipalValue1Enum); _assert_(sigma1_input);
	Input* sigma2_input = element->GetInput(DamageEffectiveStressPrincipalValue2Enum); _assert_(sigma2_input);
	Input* sigma3_input = element->GetInput(DamageEffectiveStressPrincipalValue3Enum); _assert_(sigma3_input);
	Input* tau1_input   = element->GetInput(DamageEffectiveDeviatoricStressPrincipalValue1Enum); _assert_(tau1_input);
	Input* tau2_input   = element->GetInput(DamageEffectiveDeviatoricStressPrincipalValue2Enum); _assert_(tau2_input);
	Input* tau3_input   = element->GetInput(DamageEffectiveDeviatoricStressPrincipalValue3Enum); _assert_(tau3_input);

	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		sigma1_input->GetInputValue(&s1,gauss);
		sigma2_input->GetInputValue(&s2,gauss);
		sigma3_input->GetInputValue(&s3,gauss);
		tau1_input->GetInputValue(&t1,gauss);
		tau2_input->GetInputValue(&t2,gauss);
		tau3_input->GetInputValue(&t3,gauss);

		/*Calculate stress criterion*/
		P = -(s1+s2+s3)/3; /* recover pressure */

		if(isequivstress==0){ /* von Mises */
			sigma_equiv[i]=sqrt(3/2*(t1*t1 + t2*t2 + t3*t3));
		}else if(isequivstress==1){ /* max principal stress */
			if(iseffective_stress == 0){
				sigma_equiv[i]=s1;
			}else if(iseffective_stress == 1){
				/* Eq. 16 in Keller and Hutter (2014) */
				sigma_equiv[i]=t1 - P;
			}
		}else if(isequivstress==2){ /* Hayhurst criterion */
			if(iseffective_stress == 0){
				sigma_equiv[i]=alpha*s1 + beta*sqrt(3/2*(t1*t1 + t2*t2 + t3*t3)) + (1-alpha-beta)*(s1 + s2 + s3);
			}else if(iseffective_stress== 1){
				/* Eq. 16 in Keller and Hutter (2014) */
				sigma_equiv[i]=alpha*(t1 - P) + beta*sqrt(3/2*(t1*t1 + t2*t2 + t3*t3)) + (1-alpha-beta)*(s1 + s2 + s3);
			}
		}else if(isequivstress==3){ /* Coulomb criteron */
			sigma_equiv[i]=(s1-s3)/2 + mu*(s1+s3)/2;
		}
	}

	/* Assign values */
	element->AddInput(DamageEffectiveStressEquivalentEnum,sigma_equiv,P1DGEnum);

	/* Clear memory */
	xDelete<IssmDouble>(sigma_equiv);
	delete gauss;
}/*}}}*/
