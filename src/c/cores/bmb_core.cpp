/*!\file: bmb_core.cpp
 * \brief: core of the bmb (Basal mass balance) solution 
 */ 

#include <float.h>
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
		int        stabilization;
		int        stabilizationMomentum;
		bool       ismass, ismomentum, isheat, issalt;
		IssmDouble time=0.0;
		IssmDouble timeglobal=0.0;
		IssmDouble time_year=0.0;
		IssmDouble time_day=0.0;
		IssmDouble subfinaltime; /*Global time step from ISSM*/
		IssmDouble yts;/*year to seconds*/
		IssmDouble dt; /*Time step in LADDIE*/
		int        isnonlinear;

		/*Get time step in ice and LADDIE*/
		femmodel->parameters->FindParam(&timeglobal,TimeEnum);
		femmodel->parameters->FindParam(&subfinaltime,TimesteppingTimeStepEnum);
		femmodel->parameters->FindParam(&dt,BasalforcingsLaddieSubTimestepDummyEnum);
		femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
		femmodel->parameters->FindParam(&timestepping,TimesteppingTypeEnum);
		femmodel->parameters->FindParam(&ismass,BasalforcingsLaddieIsMassEnum);
		femmodel->parameters->FindParam(&ismomentum,BasalforcingsLaddieIsMomentumEnum);
		femmodel->parameters->FindParam(&isheat,BasalforcingsLaddieIsHeatEnum);
		femmodel->parameters->FindParam(&issalt,BasalforcingsLaddieIsSaltEnum);
		femmodel->parameters->FindParam(&stabilization,BasalforcingsLaddieStabilizationMomentumEnum);

		femmodel->parameters->FindParam(&isnonlinear,BasalforcingsIsNonlinearEnum);

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

		/*Prepare ambient temperature and salinity*/
		_printf0_("Go to solve LADDIE!\n");

		if(VerboseSolution()) _printf0_("   prepare diagnositc values as zeros\n");
		LaddieDiagnosisInitialize(femmodel);

		if(VerboseSolution()) _printf0_("   preapre ambient temperature and salinity\n");
		UpdateLaddieAmbientFieldx(femmodel);

		if(VerboseSolution()) _printf0_("   prepare friction velocity\n");
		UpdateLaddieFrictionVelocityx(femmodel);

		if(VerboseSolution()) _printf0_("   prepare density and effective gravity\n");
		UpdateLaddieDensityAndEffectiveGravityx(femmodel);

		if(VerboseSolution()) _printf0_("   prepare melting rate\n");
		UpdateLaddieMeltratex(femmodel);	

		if(VerboseSolution()) _printf0_("   prepare entrainment rate\n");
		UpdateLaddieEntrainmentRatex(femmodel);

		while(time < subfinaltime - (yts*DBL_EPSILON)){
			/*Do not exceed final time of dt ISSM*/
			if(time+dt>subfinaltime - (yts*DBL_EPSILON)){
				dt=subfinaltime-time;
			}

			/*Set new step number and time step size*/
			step+=1;
			time+=dt;

			time_year = floor((time+(timeglobal-subfinaltime))/yts);
			time_day  = ((time+(timeglobal-subfinaltime))/yts - time_year)*yts/24/3600;

			femmodel->parameters->SetParam(dt,BasalforcingsLaddieSubTimestepEnum);
			_printf0_("   Laddie iteration: "<< step << "/" << ceil((subfinaltime-time)/dt)+step << \
							" time [year/days]: " << std::fixed<< setprecision(0) << time_year << " yr " << \
												 " " << std::fixed << setprecision(6) << time_day << " days\n");

			/*Step#1: Calculate mass transport model of Laddie*/
			/*Store previous step's thickness value*/
			InputDuplicatex(femmodel,BasalforcingsLaddieThicknessEnum,BasalforcingsLaddieThicknessOldEnum);
			InputDuplicatex(femmodel,BasalforcingsLaddieVxEnum,BasalforcingsLaddieVxOldEnum);
			InputDuplicatex(femmodel,BasalforcingsLaddieVyEnum,BasalforcingsLaddieVyOldEnum);

			if (ismass){
				if(VerboseSolution()) _printf0_("   computing Laddie mass transport\n");
				femmodel->SetCurrentConfiguration(BasalforcingsLaddieMassAnalysisEnum);
				solutionsequence_linear(femmodel);
			}
			/*Update thickness change (dthk/dt) depending on time*/

			/*Step#2: Calculate momentum equation of Laddie*/
			if (ismomentum){
				if(VerboseSolution()) _printf0_("   computing Laddie momentum equation\n");
				femmodel->SetCurrentConfiguration(BasalforcingsLaddieMomentumAnalysisEnum);
				if(stabilizationMomentum==4){
					solutionsequence_fct(femmodel);
				}
				else{
					//solutionsequence_newton(femmodel);
					if(isnonlinear==0) solutionsequence_nonlinear(femmodel,true);
					else if(isnonlinear==1) solutionsequence_laddie_nonlinear(femmodel);
				}
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

			/*Update all fields after laddie simulation.*/
			if(VerboseSolution()) _printf0_("   preapre ambient temperature and salinity\n");
			UpdateLaddieAmbientFieldx(femmodel);

			if(VerboseSolution()) _printf0_("   prepare friction velocity\n");
			UpdateLaddieFrictionVelocityx(femmodel);

			if(VerboseSolution()) _printf0_("   prepare density and effective gravity\n");
			UpdateLaddieDensityAndEffectiveGravityx(femmodel);

			if(VerboseSolution()) _printf0_("   prepare melting rate\n");
			UpdateLaddieMeltratex(femmodel);	

			if(VerboseSolution()) _printf0_("   prepare entrainment rate\n");
			UpdateLaddieEntrainmentRatex(femmodel);

			/*NOTE: Prepare diagnostic value for Laddie*/
			LaddieDiagnosisAdd(femmodel);
		}

		/*NOTE: Plume model contains some an osciallation, therefore, we average over the time. Divide cumulative LaddieDiagnosis variables with number of timestep*/
		LaddieDiagnosisValues(femmodel,subfinaltime);
	}

	/*Call module now*/
	FloatingiceMeltingRatex(femmodel);

	/*Extrude basal melt if not default melting rate (which may be a transient input that can't be extruded)*/
	if(basalforcing_model!=FloatingMeltRateEnum){
		femmodel->parameters->SetParam(BasalforcingsFloatingiceMeltingRateEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
	}
}
