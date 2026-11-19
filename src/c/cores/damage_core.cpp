/*
 * \brief: damage_core.cpp: core for the damage solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void damage_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(DAMAGECORE);

	/*intermediary*/
	int    damage_model;
	bool   save_results;
	bool   dakota_analysis     = false;
	int    domaintype;
	int    solution_type,stabilization;
	int    numoutputs          = 0;
	char   **requested_outputs = NULL;

	//first recover parameters common to all solutions
	femmodel->parameters->FindParam(&damage_model,DamageModelEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,DamageEvolutionNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,DamageEvolutionRequestedOutputsEnum);
	femmodel->parameters->FindParam(&stabilization,DamageStabilizationEnum);

	/*Using default damage model*/
	if (damage_model==DamageEvolutionEnum){
		DamageEvolutionAnalysis* analysis = new DamageEvolutionAnalysis();
		int damagelaw;
		femmodel->parameters->FindParam(&damagelaw,DamageLawEnum);

		/* Duplicate: */
		InputDuplicatex(femmodel,DamageDEnum,DamageDOldEnum);
		InputDuplicatex(femmodel,DamageDbarEnum,DamageDbarOldEnum);

		if(damagelaw!=0 && damagelaw!=10){
			if(VerboseSolution()) _printf0_("   update effective stress\n");
			analysis->ComputeEffectiveStress(femmodel);
			if(VerboseSolution()) _printf0_("   update stress equivalent\n");
			analysis->ComputeStressEquivalent(femmodel);
		}

		if(VerboseSolution()) _printf0_("   computing damage\n");
		Damagex(femmodel); /* optionally calculate damage analytically first */
		femmodel->SetCurrentConfiguration(DamageEvolutionAnalysisEnum);

		if(VerboseSolution()) _printf0_("   update damage input F\n");
		analysis->ComputeDamageFInput(femmodel);
		if(stabilization==3){ /* FCT */
			solutionsequence_fct(femmodel);
		}
		else{ /* Continuous Galerkin */
			solutionsequence_damage_nonlinear(femmodel);
		}
		delete analysis;

		if(domaintype==Domain2DhorizontalEnum){
			/*Duplicate DamageDbar to DamageD*/
			InputDuplicatex(femmodel,DamageDbarEnum,DamageDEnum);	
		}
		else if((domaintype==Domain3DEnum) | (domaintype==Domain2DverticalEnum)){
			/* Depth average for DamageDbar in 3D domain */
			femmodel->parameters->SetParam(DamageDEnum,InputToDepthaverageInEnum);
			femmodel->parameters->SetParam(DamageDbarEnum,InputToDepthaverageOutEnum);
			depthaverage_core(femmodel);
			
			femmodel->parameters->SetParam(DamageDbarEnum,InputToExtrudeEnum);
			extrudefrombase_core(femmodel);
		}
	}
	else{
		_error_("Damage model "<< EnumToStringx(damage_model) <<" not supported yet");
	}

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	/*Free resources:*/
	if(numoutputs){
		for(int i=0;i<numoutputs;i++){
			xDelete<char>(requested_outputs[i]);
		}
		xDelete<char*>(requested_outputs);
	}

	/*End profiler*/
	femmodel->profiler->Stop(DAMAGECORE);
}
