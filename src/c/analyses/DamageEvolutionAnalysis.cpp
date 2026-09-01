#include "./DamageEvolutionAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void DamageEvolutionAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

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
	iomodel->FetchDataToInput(inputs,elements,"md.damage.D",DamageDEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum);

	/*Initialize requested outptus in case they are not defined later for this partition*/
   iomodel->ConstantToInput(inputs,elements,0.,DamageEffectiveStressEquivalentEnum,P1Enum);
   iomodel->ConstantToInput(inputs,elements,0.,DamageEffectiveStressInvariant1Enum,P1Enum);
   iomodel->ConstantToInput(inputs,elements,0.,DamageEffectiveStressInvariant2Enum,P1Enum);
   iomodel->ConstantToInput(inputs,elements,0.,DamageEffectiveStressInvariant3Enum,P1Enum);

   iomodel->ConstantToInput(inputs,elements,0.,DamageEffectiveStressPrincipalValue1Enum,P1Enum);
   iomodel->ConstantToInput(inputs,elements,0.,DamageEffectiveStressPrincipalValue2Enum,P1Enum);
   iomodel->ConstantToInput(inputs,elements,0.,DamageEffectiveStressPrincipalValue3Enum,P1Enum);
}/*}}}*/
void DamageEvolutionAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*Intermediaries*/
	int         numoutputs;
	char**      requestedoutputs = NULL;

	/*retrieve some parameters: */
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.law",DamageLawEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.stabilization",DamageStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.max_damage",DamageMaxDamageEnum));

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.damage.requested_outputs");
	parameters->AddObject(new IntParam(DamageEvolutionNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(DamageEvolutionRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.damage.requested_outputs");

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
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.ispressure_ssa",DamageIsPressureSSAEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.isPeff",DamageIsPeffEnum));
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
void           DamageEvolutionAnalysis::PreCore(FemModel* femmodel){/*{{{*/
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
	IssmDouble tau_xx, tau_xy, tau_yy, tau_xz, tau_yz, tau_zz;
	IssmDouble s_xx,s_xy,s_xz,s_yy,s_yz,s_zz,s1,s2,stmp;
	IssmDouble s_inv1, s_inv2;
	IssmDouble pressure;
	IssmDouble Chi,Psi,PosPsi,NegPsi;
	IssmDouble damage,sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,sigma_zz,stressMaxPrincipal;
	int equivstress,domaintype,dim;

	IssmDouble  k_sigma=0.0; /* k_sigma in exponential damage model*/
	int isdamage_exponent;
	int ispressure_ssa;

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
	element->FindParam(&ispressure_ssa,DamageIsPressureSSAEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Get problem dimension*/
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
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
	Input* stress_inv1_input = element->GetInput(DamageEffectiveStressInvariant1Enum); _assert_(stress_inv1_input);

	/*retrieve the desired type of equivalent stress*/
	element->FindParam(&equivstress,DamageEquivStressEnum);

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		/* Get first invariant of stress tensor */
		damage_input->GetInputValue(&damage,gauss);

		/* Call effective stress criterion from "DamageStressEquivalentEnum" */
		stress_equivalent_input->GetInputValue(&Chi,gauss);
		stress_inv1_input->GetInputValue(&s_inv1,gauss);

		/*Calculate k_simga value*/
		if(isdamage_exponent==0){
			k_sigma = c3;
		}else if(isdamage_exponent==1){
			/* Eq 23. Pralong et al. (2005) */
			IssmDouble k1=3.75e-3; // [Pa-1] Table 3 in Pralong et al. (2005)

			k_sigma = k1 * ( sqrt(max(s_inv1,0.0)) - healing * sqrt(max(-s_inv1,0.0)) );
		}else if(isdamage_exponent==2){
			/*FIXME: Hard coding. Parameters from Duddu et al. (2020)*/
			IssmDouble k1=-2.63; // unit: -
			IssmDouble k2=7.24*1e-6; // unit:  Pa-1

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
	IssmDouble  vel,vx,vy,vz,dvxdx,dvydy,dvzdz,dvx[3],dvy[3],dvz[3];
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
	if(dim==3){
		vz_input=element->GetInput(VzEnum); _assert_(vz_input);
	}

	if(dim==2) h=element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

		if(dim==3){
			vz_input->GetInputValue(&vz,gauss);
			vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
		}

		/*Transient term*/
		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j];

		dvxdx=dvx[0];
		dvydy=dvy[1];
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
			dvzdz=dvz[2];
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*\phi_i \phi_j \nabla\cdot v*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy+dvzdz);
					/*\phi_i v\cdot\nabla\phi_j*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j] + vz*dbasis[2*numnodes+j]);
				}
			}
		}

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
			else{
				/*Streamline upwinding*/
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				D[0*dim+0]=h/(2.0*vel)*vx*vx;
				D[1*dim+0]=h/(2.0*vel)*vy*vx;
				D[0*dim+1]=h/(2.0*vel)*vx*vy;
				D[1*dim+1]=h/(2.0*vel)*vy*vy;
			}
		}
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

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		newdamage[i]=solution[doflist[i]];
		/*Check solution*/
		if(xIsNan<IssmDouble>(newdamage[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(newdamage[i])) _error_("Inf found in solution vector");
		/*Enforce D < max_damage and D > 0 */
		if(newdamage[i]>max_damage) newdamage[i]=max_damage;
		else if(newdamage[i]<0.)    newdamage[i]=0.;
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
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented yet");
	}
	
	/* Intermediaries */
	IssmDouble tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz; /* Deviatoric stress tensor components [Pa] */
	IssmDouble sigma_xx, sigma_xy, sigma_xz, sigma_yy, sigma_yz, sigma_zz; /* Cauchy stress tensor components [Pa] */
	IssmDouble s1, s2, s3; /* principal effective stress [Pa] */
	IssmDouble stmp;
	IssmDouble D; /* Damage value [unit: 1]*/
	IssmDouble z; /* elevation [m] */

	IssmDouble rho_w, g;
	IssmDouble *sigma_equiv = NULL; /* equivalent stress [Pa] */
	IssmDouble *sigma_inv1 = NULL; /* first invariant stress [Pa] */
	IssmDouble *sigma_inv2 = NULL; /* second invariant stress [Pa] */
	IssmDouble *sigma_inv3 = NULL; /* third invariant stress [Pa] */
	IssmDouble alpha, beta; /* parameter for Hayhurst */
	IssmDouble mu; /* parameter for Coulomb */

	/*    Principal stress for each vertex */
	IssmDouble *sigma_1 = NULL;
	IssmDouble *sigma_2 = NULL;
	IssmDouble *sigma_3 = NULL;

	IssmDouble *xyz_list = NULL;
	
	/*    ice flow equation */
	bool isSSA, isL1L2, isMOLHO, isHO, isFS;
	
	/* Pressure term
	P  : pressure for stress tensor [Pa]
	Pi : Overburden pressure [Pa] 
	Pw : Water pressure [Pa]
	*/
	IssmDouble P, Pi, Pw;
	bool isPeff; /* if true, we compute the effective stress, else we compute the Cauchy stress */

	int isequivstress; /* stress equivalent stress */
	int ispressure_ssa; /* treatment for pressure in SSA2D model */

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	sigma_equiv = xNewZeroInit<IssmDouble>(numnodes);

	sigma_inv1  = xNewZeroInit<IssmDouble>(numnodes);
	sigma_inv2  = xNewZeroInit<IssmDouble>(numnodes);
	sigma_inv3  = xNewZeroInit<IssmDouble>(numnodes);

	sigma_1  = xNewZeroInit<IssmDouble>(numnodes);
	sigma_2  = xNewZeroInit<IssmDouble>(numnodes);
	sigma_3  = xNewZeroInit<IssmDouble>(numnodes);

	/* Retrieve constants */
	element->FindParam(&rho_w,MaterialsRhoSeawaterEnum);
	element->FindParam(&g,ConstantsGEnum);
	element->FindParam(&isequivstress,DamageEquivStressEnum);
	element->FindParam(&alpha,DamageEquivStressAlphaEnum);
	element->FindParam(&beta,DamageEquivStressBetaEnum);
	element->FindParam(&mu,DamageEquivStressMuEnum);
	element->FindParam(&ispressure_ssa,DamageIsPressureSSAEnum);
	element->FindParam(&isPeff,DamageIsPeffEnum);

	element->FindParam(&isSSA,FlowequationIsSSAEnum);
	element->FindParam(&isL1L2,FlowequationIsL1L2Enum);
	element->FindParam(&isMOLHO,FlowequationIsMOLHOEnum);
	element->FindParam(&isHO,FlowequationIsHOEnum);
	element->FindParam(&isFS,FlowequationIsFSEnum);

	/* Retrieve all input and parameters */
	element->GetVerticesCoordinates(&xyz_list);

	/* Precompute deviatoric stress tensor*/
	/* NOTE: ComputeDeviatoricStressTensor already contains damage i.g. tau_eff = tau/(1-D)*E; therefore, tau_eff is the effective deviatoric stress. */
	element->ComputeDeviatoricStressTensor();
	// if(dim==3){
	// 	/*Only works in 3d because the pressure is defined*/
	// 	element->StressMaxPrincipalCreateInput();
	// }

	/* Retrieve what we need: */
	Input* tau_xx_input  = element->GetInput(DeviatoricStressxxEnum);     _assert_(tau_xx_input);
	Input* tau_xy_input  = element->GetInput(DeviatoricStressxyEnum);     _assert_(tau_xy_input);
	Input* tau_yy_input  = element->GetInput(DeviatoricStressyyEnum);     _assert_(tau_yy_input);
	Input* tau_xz_input  = NULL;
	Input* tau_yz_input  = NULL;
	Input* tau_zz_input  = NULL;
	if(dim==3){
		tau_xz_input  = element->GetInput(DeviatoricStressxzEnum);     _assert_(tau_xz_input);
		tau_yz_input  = element->GetInput(DeviatoricStressyzEnum);     _assert_(tau_yz_input);
		tau_zz_input  = element->GetInput(DeviatoricStresszzEnum);     _assert_(tau_zz_input);
	}
	/* NOTE: pressure in model would be overburden pressure */
	Input* pressure_input = element->GetInput(PressureEnum);          _assert_(pressure_input);

	Input* damage_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
	}

	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);

		damage_input->GetInputValue(&D,gauss);

		tau_xx_input->GetInputValue(&tau_xx,gauss);
		tau_xy_input->GetInputValue(&tau_xy,gauss);
		tau_yy_input->GetInputValue(&tau_yy,gauss);
		if(dim==3){
			tau_xz_input->GetInputValue(&tau_xz,gauss);
			tau_yz_input->GetInputValue(&tau_yz,gauss);
			tau_zz_input->GetInputValue(&tau_zz,gauss);
			// stressMaxPrincipal_input = element->GetInput(StressMaxPrincipalEnum); _assert_(stressMaxPrincipal_input);
		}
		
		/* Compute pressure in Cauchy stress tensor*/
		pressure_input->GetInputValue(&Pi,gauss);
		/* Eq. 15 in Huth et al. (2021) / SSA approximation with hydrostatic assumption (Greve and Blatter et al., 2009) */
		if (isSSA || isHO || isMOLHO) {
			/* NOTE: deviratoric stress should be damaged? ~ tau^{D} = tau/(1-D) */
			Pi = Pi - tau_xx - tau_yy; 
		}else{
			_error_("not implemented yet");
		}
	
		if(dim==2){
			P = 0.0; /* pressure at surface*/
			if(ispressure_ssa==0){ /* Assume pressure at surface */
				P=0.0;
			}else if(ispressure_ssa==1){
				P=0.5*Pi; /* P is defined in 2d as half of the 3d P */
			}else if(ispressure_ssa==2){
				/*Nothing to be done*/
				P=1.0*Pi; /* pressure at bed */
			}
		}else if(dim==3){
			z = xyz_list[i*3+2];
			Pw = 0.0;
			if (isPeff){
				if (z > 0){ /* explicitly set Pw value */
					Pw = 0.0;
				}else if(z <= 0){
					Pw = rho_w*g*(0-z);
				}
			}
			P = Pi - Pw;
		}

		/* Compute effective Cauchy stress tensor baed on deviatoric stress */
		//sigma_xx =   tau_xx/(1-D) - P;
		//sigma_xy =   tau_xy/(1-D);
		//sigma_yy =   tau_yy/(1-D) - P;
		//if(dim==3){
		//	sigma_xz = tau_xz/(1-D);
		//	sigma_yz = tau_yz/(1-D);
		//	sigma_zz = tau_zz/(1-D) - P;
		//}
		sigma_xx =   tau_xx - P;
		sigma_xy =   tau_xy;
		sigma_yy =   tau_yy - P;
		if(dim==3){
			sigma_xz = tau_xz;
			sigma_yz = tau_yz;
			sigma_zz = tau_zz - P;
		}

		/*Calculate principal effective stresses*/
		if(dim==2){
			/* Compute principal effective stresses */
			s3 = -P;
			Matrix2x2Eigen(&s1,&s2,NULL,NULL,sigma_xx,sigma_xy,sigma_yy);

			if(s2>s1){stmp=s2; s2=s1; s1=stmp;}

			if(isequivstress==0){ /* von Mises */
				sigma_equiv[i]=sqrt(s1*s1-s1*s2+s2*s2);
			}else if(isequivstress==1){ /* max principal stress */
				sigma_equiv[i]=s1;
			}
			else if(isequivstress==2){ /* Hayhurst stress invariant */
				sigma_equiv[i]=alpha*s1 + beta*sqrt(s1*s1-s1*s2+s2*s2) + (1-alpha-beta)*(s1 + s2);
			}else if(isequivstress==3){ /* Coulomb criteron */
				sigma_equiv[i]=(s1-s3)/2 + mu*(s1+s2)/2;
			}
		}else if(dim==3){
			/* Compute principal effective stresses*/
			/* FIXME: Now, stress equivalent would be only computed for SSA and HO, not FS.*/
			if(isSSA || isHO){ /* Hydrostatic assumption : neglect sigma_zx, sigma_zy */
				Matrix3x3Eigen(&s1,&s2,&s3,NULL,NULL,NULL,
					sigma_xx,sigma_xy,sigma_xz,
					sigma_xy,sigma_yy,sigma_yz,
					0.0, 0.0, sigma_zz);
			}else{
				_error_("not implemented yet");
			}

			/* Ordering large value (descending order) ... */
			if(s1<s2); swap(s1,s2);
			if(s1<s3); swap(s1,s3);
			if(s2<s3); swap(s2,s3);
		
			//inv1 = sigma_xx + sigma_yy + sigma_zz;
			//inv2 = sigma_xx*sigma_yy + sigma_yy*sigma_zz + sigma_zz*sigma_xx - sigma_xy*sigma_xy - sigma_yz*sigma_yz - sigma_xz*sigma_xz;

			if(isequivstress==0){ /* von Mises */
				sigma_equiv[i]=sqrt(((s1-s2)*(s1-s2)+(s2-s3)*(s2-s3)+(s3-s1)*(s3-s1))/2.);
			}else if(isequivstress==1){ /* max principal stress */
				sigma_equiv[i]=s1;
			}else if(isequivstress==2){ /* Hayhurst stress invariant */
				sigma_equiv[i]=alpha*s1 + beta*sqrt(((s1-s2)*(s1-s2)+(s2-s3)*(s2-s3)+(s3-s1)*(s3-s1))/2.) + (1-alpha-beta)*(s1 + s2 + s3);
			}else if(isequivstress==3){ /* Coulomb criteron */
				sigma_equiv[i]=(s1-s3)/2 + mu*(s1+s3)/2;
			}
		}
		
		/* Now assign invariant values from principal stresses */
		sigma_inv1[i] = s1 + s2 + s3;
		sigma_inv2[i] = s1*s2 + s1*s3 + s2*s3;
		sigma_inv3[i] = s1*s2*s3;

		sigma_1[i] = s1;
		sigma_2[i] = s2;
		sigma_3[i] = s3;
	}

	/* Assign values */
	element->AddInput(DamageEffectiveStressEquivalentEnum,sigma_equiv,P1Enum);
	element->AddInput(DamageEffectiveStressInvariant1Enum,sigma_inv1,P1Enum);
	element->AddInput(DamageEffectiveStressInvariant2Enum,sigma_inv2,P1Enum);
	element->AddInput(DamageEffectiveStressInvariant3Enum,sigma_inv3,P1Enum);

	element->AddInput(DamageEffectiveStressPrincipalValue1Enum,sigma_1,P1Enum);
	element->AddInput(DamageEffectiveStressPrincipalValue2Enum,sigma_2,P1Enum);
	element->AddInput(DamageEffectiveStressPrincipalValue3Enum,sigma_3,P1Enum);

	/* Clear memory */
	xDelete<IssmDouble>(sigma_equiv);
	xDelete<IssmDouble>(sigma_inv1);
	xDelete<IssmDouble>(sigma_inv2);
	xDelete<IssmDouble>(sigma_inv3);
	xDelete<IssmDouble>(sigma_1);
	xDelete<IssmDouble>(sigma_2);
	xDelete<IssmDouble>(sigma_3);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
}/*}}}*/
