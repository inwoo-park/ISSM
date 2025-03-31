/*!\file: bmb_core.cpp
 * \brief: core of the bmb (Basal mass balance) solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void bmb_core(FemModel* femmodel){

	/*First, get BMB model from parameters*/
	int  basalforcing_model;
	bool isplume = false;
	femmodel->parameters->FindParam(&basalforcing_model,BasalforcingsEnum);

	if(VerboseSolution()) _printf0_("   computing basal mass balance\n");

	/*In some cases we need to run additional analyses to get the required input data*/
	if(basalforcing_model==BasalforcingsPicoEnum){
		femmodel->parameters->FindParam(&isplume,BasalforcingsPicoIsplumeEnum);
		if(isplume){
			femmodel->SetCurrentConfiguration(L2ProjectionBaseAnalysisEnum);
			femmodel->parameters->SetParam(BaseSlopeXEnum,InputToL2ProjectEnum);
			solutionsequence_linear(femmodel);
			femmodel->parameters->SetParam(BaseSlopeYEnum,InputToL2ProjectEnum);
			solutionsequence_linear(femmodel);
			femmodel->SetCurrentConfiguration(GLheightadvectionAnalysisEnum);
			solutionsequence_linear(femmodel);
		}
	}
	else if(basalforcing_model==BasalforcingsLaddieEnum){
		/*Sub-ice shelf melting with LADDIE simulation*/
		int        timestepping; /*check TimeSteppingEnum*/
		int        step=0;
		bool       ismass, ismomentum, isheat, issalt;
		IssmDouble time=0.0;
		IssmDouble subfinaltime; /*Global time step from ISSM*/
		IssmDouble yts;/*year to seconds*/
		IssmDouble dt; /*Time step in LADDIE*/

		/*Get time step in ice and LADDIE*/
		femmodel->parameters->FindParam(&subfinaltime,TimesteppingTimeStepEnum);
		femmodel->parameters->FindParam(&dt,BasalforcingsLaddieSubTimestepDummyEnum);
		femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
		femmodel->parameters->FindParam(&timestepping,TimesteppingTypeEnum);
		femmodel->parameters->FindParam(&ismass,BasalforcingsLaddieIsMassEnum);
		femmodel->parameters->FindParam(&ismomentum,BasalforcingsLaddieIsMomentumEnum);
		femmodel->parameters->FindParam(&isheat,BasalforcingsLaddieIsHeatEnum);
		femmodel->parameters->FindParam(&issalt,BasalforcingsLaddieIsSaltEnum);

		/*Check time stepping type*/
		switch(timestepping){
			case FixedTimesteppingEnum:
				/*Nothing to do*/
				break;
			case AdaptiveTimesteppingEnum:
				_error_("Time stepping \""<<EnumToStringx(timestepping)<<"\" not supported yet");
				break;
			default:
				_error_("Time stepping \""<<EnumToStringx(timestepping)<<"\" not supported yet");
		}

		/*Step#1: prepare ambient temperature and salinity*/
		if(VerboseSolution()) _printf0_("      step1: Preapre ambient temperature and salinity\n");
		UpdateLaddieAmbientFieldx(femmodel);

		/*Step#2: update friction velocity*/
		if(VerboseSolution()) _printf0_("      step2: Prepare frictinon velocity\n");
		UpdateLaddieFrictionVelocityx(femmodel);

		/*Update density and effective gravity: */
		if(VerboseSolution()) _printf0_("      step3: Prepare density and effective gravity\n");
		UpdateLaddieDensityAndEffectiveGravityx(femmodel);

		/*First, initialize guess sub-ice shelf melting and entrainment rate.*/
		if(VerboseSolution()) _printf0_("   computing melting rate\n");
		FloatingiceMeltingRateLaddiex(femmodel);	

		if(VerboseSolution()) _printf0_("   computing entrainment rate\n");
		UpdateLaddieEntrainmentRatex(femmodel);

		_printf0_("Go to solve LADDIE!\n");
		while(time < subfinaltime - (yts*DBL_EPSILON)){
			/*Do not exceed final time of dt ISSM*/
			if(time+dt>subfinaltime - (yts*DBL_EPSILON)){
				dt=subfinaltime-time;
			}

			/*Set new step number and time step size*/
			step+=1;
			time+=dt;
			femmodel->parameters->SetParam(dt,BasalforcingsLaddieSubTimestepEnum);
			_printf0_("   Laddie iteration: "<< step << "/" << ceil((subfinaltime-time)/dt)+step << \
							" time [days]: " << time/24/3600 << "\n");
			//if(VerboseSolution()) _printf0_("   Laddie time: "<<time/24/3600<<" days\n");

			/*Step#1: Calculate mass transport model of Laddie*/
			if (ismass){
				if(VerboseSolution()) _printf0_("   computing Laddie mass transport\n");
				femmodel->SetCurrentConfiguration(BasalforcingsLaddieMassAnalysisEnum);
				solutionsequence_linear(femmodel);
			}

			/*Step#2: Calculate momentum equation of Laddie*/
			if (ismomentum){
				if(VerboseSolution()) _printf0_("   computing Laddie momentum equation\n");
				femmodel->SetCurrentConfiguration(BasalforcingsLaddieMomentumAnalysisEnum);
				solutionsequence_linear(femmodel);
			}

			/*Step#3: Calculate heat equation of Laddie*/
			if (isheat){
				if(VerboseSolution()) _printf0_("   computing Laddie heat equation\n");
				femmodel->SetCurrentConfiguration(BasalforcingsLaddieHeatAnalysisEnum);
				solutionsequence_linear(femmodel);
			}

			/*Step#4: Calculate salt equation of Laddie*/
			if (issalt){
				if(VerboseSolution()) _printf0_("   computing Laddie salt equation\n");
				femmodel->SetCurrentConfiguration(BasalforcingsLaddieSaltAnalysisEnum);
				solutionsequence_linear(femmodel);
			}

			/*Update friction velocity*/
			UpdateLaddieFrictionVelocityx(femmodel);

			/*Update density and effective gravity*/
			UpdateLaddieDensityAndEffectiveGravityx(femmodel);

			/*First, initialize guess sub-ice shelf melting and entrainment rate.*/
			if(VerboseSolution()) _printf0_("   computing melting rate and entrainment rate\n");
			FloatingiceMeltingRateLaddiex(femmodel);	
			UpdateLaddieEntrainmentRatex(femmodel);
		}
	}

	/*Call module now*/
	FloatingiceMeltingRatex(femmodel);

	/*Extrude basal melt if not default melting rate (which may be a transient input that can't be extruded)*/
	if(basalforcing_model!=FloatingMeltRateEnum){
		femmodel->parameters->SetParam(BasalforcingsFloatingiceMeltingRateEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
	}
}
