/*!\file FloatingiceMeltingRateLaddiex
 * \brief: calculates Floating ice melting rate
 */

#include "./FloatingiceMeltingRateLaddiex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "./../../classes/Inputs/DatasetInput.h"
#include "../InputDuplicatex/InputDuplicatex.h"

void FloatingiceMeltingRateLaddiex(FemModel* femmodel){/*{{{*/
	/*Duplicate time-averaged floating ice melting rate as BasalforcingsFloatingiceMeltingRateEnum*/
	
	//InputDuplicatex(femmodel,BasalforcingsLaddieMeltingRateEnum,BasalforcingsFloatingiceMeltingRateEnum);
	InputDuplicatex(femmodel,BasalforcingsLaddieDiagMeltingRateEnum,BasalforcingsFloatingiceMeltingRateEnum);
}/*}}}*/

IssmDouble GetDensityDifferencex(IssmDouble rho0, IssmDouble T, IssmDouble S, IssmDouble Ta, IssmDouble Sa){/*{{{*/
	/*
	 Explain
	 Calculate the density difference between the mixed layer and ambient ocean below the layer.
	 */
	IssmDouble drho;
	IssmDouble alpha=3.733e-5; /*thermal expansion coefficient [degC-1]*/
   IssmDouble beta=7.843e-4; /*haline contraction coefficient [psu-1]*/

	drho = rho0*(-alpha*(Ta-T) + beta*(Sa-S));

	return drho;
}/*}}}*/
IssmDouble GetEffectiveGravitationAccelerationx(IssmDouble g, IssmDouble rho0, IssmDouble T, IssmDouble S, IssmDouble Ta, IssmDouble Sa){/*{{{*/
    /*
    Explain
    Calculate the effective gravitational density due to density difference in seawater.

    Inputs
    * rho0: initial ocean density [kg m-3]
    * Ta: ambient ocean temperature [degC]
    * Sa: ambient ocean salinity [psu]
    * T: ocean temprature in mixed layer [degC]
    * S: ocean salinity in mixed layer [psu]
    
    Outputs
    * g_e: effective gravitational acceleration (m s-2)

    See also Eq. 6 and Eq. 7 in Lambert et al. (2023).
    */
	 //IssmDouble g; /*Gravitational acceleration*/
    IssmDouble g_e; /*Effective gravitational acceleration*/
    IssmDouble drho; /*Density difference between the layer and the ambient water below the layer*/
    IssmDouble alpha=3.733e-5; /*thermal expansion coefficient [degC-1]*/
    IssmDouble beta=7.843e-4; /*haline contraction coefficient [psu-1]*/

    /*Calculate densitify difference*/
    drho = rho0*(-alpha*(Ta-T) + beta*(Sa-S));

    /*Calculate effective gravitational density*/
    g_e = g*drho/rho0;

    return g_e;
}/*}}}*/

void LaddieDiagnosisInitialize(FemModel *femmodel){/*{{{*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		IssmDouble* zeros=xNew<IssmDouble>(numvertices);

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			zeros[iv] = 0.0;
		}

		/*Assign values*/
		element->AddInput(BasalforcingsLaddieDiagMeltingRateEnum,zeros,P1DGEnum);
		//element->AddInput(BasalforcingsLaddieDiagVxEnum,zeros,P1DGEnum);
		//element->AddInput(BasalforcingsLaddieDiagVyEnum,zeros,P1DGEnum);
		//element->AddInput(BasalforcingsLaddieDiagVelEnum,zeros,P1DGEnum);
		//element->AddInput(BasalforcingsLaddieDiagThicknessEnum,zeros,P1DGEnum);
		//element->AddInput(BasalforcingsLaddieDiagTemperatureEnum,zeros,P1DGEnum);
		//element->AddInput(BasalforcingsLaddieDiagSalinityEnum,zeros,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(zeros);
		delete gauss;
	}
}/*}}}*/
void LaddieDiagnosisAdd(FemModel *femmodel){/*{{{*/
	/*
	Accumulate calculated sub-ice shelf melting rate for smooth-out oscillation in plume model.
	 */

	IssmDouble subtimestep;

	/*Retrieve all inputs*/
	femmodel->parameters->FindParam(&subtimestep,BasalforcingsLaddieSubTimestepEnum);

	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		Input* melt_input     = element->GetInput(BasalforcingsLaddieMeltingRateEnum); _assert_(melt_input);
		Input* melt_diag_input= element->GetInput(BasalforcingsLaddieDiagMeltingRateEnum); _assert_(melt_diag_input);

		IssmDouble* melt=xNew<IssmDouble>(numvertices);
		IssmDouble* melt_diag=xNew<IssmDouble>(numvertices);

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			/*Get all dataset*/
			melt_input->GetInputValue(&melt[iv],gauss);
			melt_diag_input->GetInputValue(&melt_diag[iv],gauss);

			/*Add current method*/
			melt_diag[iv] += melt[iv]*subtimestep;
		}

		/*Assign values*/
		element->AddInput(BasalforcingsLaddieDiagMeltingRateEnum,melt_diag,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(melt);
		xDelete<IssmDouble>(melt_diag);
		delete gauss;
	}
}/*}}}*/
void LaddieDiagnosisValues(FemModel *femmodel, IssmDouble timestep){/*{{{*/
	/*
		Calculate mean values during sub-time stepping in Laddie simulation.
	 */

	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		Input* melt_diag_input= element->GetInput(BasalforcingsLaddieDiagMeltingRateEnum); _assert_(melt_diag_input);

		IssmDouble* melt_diag=xNew<IssmDouble>(numvertices);

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			/*Get all dataset*/
			melt_diag_input->GetInputValue(&melt_diag[iv],gauss);

			/*Add current method*/
			melt_diag[iv] = melt_diag[iv]/timestep;
		}

		/*Assign values*/
		element->AddInput(BasalforcingsLaddieDiagMeltingRateEnum,melt_diag,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(melt_diag);
		delete gauss;
	}
}/*}}}*/

void UpdateLaddieDThicknessDtx(FemModel* femmodel){/*{{{*/
	/*
	Prepare (d thickness / dt ) for heat and salt analysis
	 */

	IssmDouble  dt;
	IssmDouble  thk1, thk0;

	IssmDouble *dthkdt;
	IssmDouble *values;

	Gauss      *gauss;

	/*Retrieve all inputs: */
	femmodel->parameters->FindParam(&dt,BasalforcingsLaddieSubTimestepEnum);

	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieDThicknessDtEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		Input *thk1_input = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thk1_input);
		Input *thk0_input = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thk0_input);
		dthkdt = xNew<IssmDouble>(numvertices);

		gauss=element->NewGauss();
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			thk1_input->GetInputValue(&thk1,gauss);
			thk0_input->GetInputValue(&thk0,gauss);
			
			/*calculate dthk/dt value: */
			dthkdt[iv] = (thk1-thk0)/dt;
		}

		/*Assign value in: */
		element->AddInput(BasalforcingsLaddieDThicknessDtEnum,dthkdt,P1DGEnum);

		xDelete<IssmDouble>(dthkdt);
		delete gauss;
	}
}/*}}}*/
void UpdateLaddieAmbientFieldx(FemModel* femmodel){/*{{{*/
	/*
	Make 3D ambient ocean temperature/salinity to 2D ambient ocean temperature/salinity considering ice base elevation and plume thickness.

	* zb: base elevation of ice shelf.
	* D: plume thickness.

	See also
	laddie/src/physics.py > "def update_ambientfields"
	 */
	
	IssmDouble *values;
	IssmDouble *forcing_depth=NULL;
	int num_depths;

	femmodel->parameters->FindParam(&forcing_depth,&num_depths,BasalforcingsLaddieForcingDepthEnum); _assert_(forcing_depth);

	/*Get ambient ocean temperature and salinity at each ice shelf point - linearly interpolate in depth and time*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		//if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
		//	values = xNewZeroInit<IssmDouble>(numvertices);
		//	element->AddInput(BasalforcingsLaddieAmbientTemperatureEnum,values,P1DGEnum);
		//	xDelete<IssmDouble>(values);

		//	values = xNewZeroInit<IssmDouble>(numvertices);
		//	element->AddInput(BasalforcingsLaddieAmbientSalinityEnum,values,P1DGEnum);
		//	xDelete<IssmDouble>(values);
		//	continue;
		//}
		
		/*Get ambient ocean temeprature and salinity on all vertices*/
		IssmDouble *Ttmp           = xNew<IssmDouble>(numvertices);
		IssmDouble *Stmp           = xNew<IssmDouble>(numvertices);
		IssmDouble *depth_vertices = xNew<IssmDouble>(numvertices);
		IssmDouble *thickness_vertices = xNew<IssmDouble>(numvertices);
		/*
		 * tf: temperature forcing
		 * sf: salinity forcing
		 */
		DatasetInput* tf_input = element->GetDatasetInput(BasalforcingsLaddieForcingTemperatureEnum); _assert_(tf_input);
		DatasetInput* sf_input = element->GetDatasetInput(BasalforcingsLaddieForcingSalinityEnum); _assert_(sf_input);

		element->GetInputListOnVertices(&depth_vertices[0],BaseEnum);
		// NOTE: If the vertices are located in grounded ice, thickness_vertices are assigned zeros, as plume thickness over grounded ice is set to zero.
		element->GetInputListOnVertices(&thickness_vertices[0],BasalforcingsLaddieThicknessEnum);

		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			IssmDouble depth=(depth_vertices[iv]-thickness_vertices[iv]);
			int offset;
			int found=binary_search(&offset,depth,forcing_depth,num_depths);
			if(!found) _error_("depth not found");

			if (offset==-1){
				/*get values for the first depth: */
				_assert_(depth<=forcing_depth[0]);
				tf_input->GetInputValue(&Ttmp[iv],gauss,0);
				sf_input->GetInputValue(&Stmp[iv],gauss,0);
			}
			else if(offset==num_depths-1){
				/*get values for the last time: */
				_assert_(depth>=forcing_depth[num_depths-1]);
				tf_input->GetInputValue(&Ttmp[iv],gauss,num_depths-1);
				sf_input->GetInputValue(&Stmp[iv],gauss,num_depths-1);
			}
			else {
				/*get values between two times [offset:offset+1], Interpolate linearly*/
				_assert_(depth>=forcing_depth[offset] && depth<forcing_depth[offset+1]);
				IssmDouble deltaz=forcing_depth[offset+1]-forcing_depth[offset];
				IssmDouble alpha2=(depth-forcing_depth[offset])/deltaz;
				IssmDouble alpha1=(1.-alpha2);
				IssmDouble tf1,tf2;
				IssmDouble sf1,sf2;
				tf_input->GetInputValue(&tf1,gauss,offset);
				tf_input->GetInputValue(&tf2,gauss,offset+1);
				Ttmp[iv] = alpha1*tf1 + alpha2*tf2;

				sf_input->GetInputValue(&sf1,gauss,offset);
				sf_input->GetInputValue(&sf2,gauss,offset+1);
				Stmp[iv] = alpha1*sf1 + alpha2*sf2;
			}
		}

		element->AddInput(BasalforcingsLaddieAmbientTemperatureEnum,Ttmp,P1DGEnum);
		element->AddInput(BasalforcingsLaddieAmbientSalinityEnum,Stmp,P1DGEnum);

		/*Clear memory: */
		xDelete<IssmDouble>(Ttmp);
		xDelete<IssmDouble>(Stmp);
		xDelete<IssmDouble>(depth_vertices);
		xDelete<IssmDouble>(thickness_vertices);
		delete gauss;
	}

	/*Clear memory: */
	xDelete<IssmDouble>(forcing_depth);
}/*}}}*/
void UpdateLaddieDensityAndEffectiveGravityx(FemModel* femmodel){/*{{{*/
	/*
	Update ocean density and effective gravitational acceleration due to buoyancy.
	
	See also
	Eq. 6 and Eq. 7 in Lambert et al. (2023): https://tc.copernicus.org/articles/17/3203/2023/
	 */

	IssmDouble  rho0; /*default ocean density [kg m-3]*/
	IssmDouble  g; /*gravitaional acceleration [m s-1]*/
	IssmDouble  T,S; /*plume temperature and salinity*/
	IssmDouble  Ta, Sa; /*ambient ocean forcing for temperature [degC] and salinity [psu]*/
	IssmDouble  alpha=3.733e-5; /*thermal expansion coefficient [degC -1]*/
	IssmDouble  beta=7.843e-4; /*haline contraction coefficient [psu-1]*/
	IssmDouble  mindrho=0.005; /*minmum density difference [kg m-3]*/
	int         convop;

	IssmDouble *values;
	int         iscatch=0;

	femmodel->parameters->FindParam(&convop,BasalforcingsLaddieConvOptionEnum);
	femmodel->parameters->FindParam(&rho0,MaterialsRhoSeawaterEnum);
	femmodel->parameters->FindParam(&g,ConstantsGEnum);

	/*Get basal friction coefficient*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set DRho and ga values to 0 if non-floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			values = xNewZeroInit<IssmDouble>(numvertices);
			for(int i=0; i<numvertices;i++) values[i] = mindrho;
			element->AddInput(BasalforcingsLaddieDRhoEnum,&values[0],P1DGEnum);
			xDelete<IssmDouble>(values);

			values = xNewZeroInit<IssmDouble>(numvertices);
			for(int i=0; i<numvertices;i++) values[i] = g*mindrho/rho0;
			element->AddInput(BasalforcingsLaddieAmbientGEnum,&values[0],P1DGEnum);
			xDelete<IssmDouble>(values);

			continue;
		}

		/*Get inputs*/
		/*Plume temperature and salinity*/
		Input* T_input =element->GetInput(BasalforcingsLaddieTEnum); _assert_(T_input);
		Input* S_input =element->GetInput(BasalforcingsLaddieSEnum); _assert_(S_input);
		/*Ambient ocean temperature and salinity*/
		Input* Ta_input=element->GetInput(BasalforcingsLaddieAmbientTemperatureEnum); _assert_(Ta_input);
		Input* Sa_input=element->GetInput(BasalforcingsLaddieAmbientSalinityEnum); _assert_(Sa_input);

		/*Initialize variable*/
		IssmDouble* ga  =xNew<IssmDouble>(numvertices);
		IssmDouble* drho=xNew<IssmDouble>(numvertices);
		IssmDouble* Tnew=xNew<IssmDouble>(numvertices);
		IssmDouble* Snew=xNew<IssmDouble>(numvertices);


		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			
			/*Prepare value*/
			T_input->GetInputValue(&T, gauss);
			S_input->GetInputValue(&S, gauss);
			Ta_input->GetInputValue(&Ta, gauss);
			Sa_input->GetInputValue(&Sa, gauss);

			/*Calculate density difference*/
			Tnew[iv] = T;
			Snew[iv] = S;
			drho[iv] = rho0*(-alpha*(Ta-T) + beta*(Sa-S));
			if (convop == 0){
				/*Prescribe minimum stratification*/
				drho[iv] = max(drho[iv], mindrho);
			}else if (convop == 1){
				/*Apply instaneous convection*/
				if (drho[iv] < mindrho){
					Tnew[iv] = Ta;
					Snew[iv] = Sa-mindrho/(rho0*beta);
					drho[iv] = rho0*(-alpha*(Ta-Tnew[iv]) + beta*(Sa-Snew[iv]));
				}
			}else{
				_error_("convOption = " << convop << " is not supported.\n");
			}
			ga[iv] = g*drho[iv]/rho0;
		}

		/*Assign value in: */
		element->AddInput(BasalforcingsLaddieDRhoEnum,&drho[0],P1DGEnum);
		element->AddInput(BasalforcingsLaddieAmbientGEnum,&ga[0],P1DGEnum);

		element->AddInput(BasalforcingsLaddieTEnum,&Tnew[0],P1Enum);
		element->AddInput(BasalforcingsLaddieSEnum,&Snew[0],P1Enum);
		
		/*Clear memory*/
		xDelete<IssmDouble>(ga);
		xDelete<IssmDouble>(drho);
		xDelete<IssmDouble>(Tnew);
		xDelete<IssmDouble>(Snew);
	}
}/*}}}*/

void UpdateLaddieFrictionVelocityx(FemModel* femmodel){/*{{{*/
	/*
	 Calculate the friction velocity Ustar defined in Jenkins et al. (2010)

	 Ustar = sqrt(Cd_top (vx^2 + vy^2 + Utide^2))

	 Inputs
	 * Cd_top: drag coefficient.
	 * vx, vy: plume velocity [m s-1].
	 * Utide: time-mean tidal velocity [m s-1].

	 Outputs
	 * Ustar - element friction velocity

	 See also
	 Eq. 13 in Lambert et al. (2023): https://tc.copernicus.org/articles/17/3203/2023/
	 */
	IssmDouble Cdtop; /*Friction coefficient [m s-1]*/
	IssmDouble vx, vy; /*plume velocity vector [m s-1]*/
	IssmDouble utide; /*tide velocity*/

	femmodel->parameters->FindParam(&Cdtop, BasalforcingsLaddieCdTopEnum);
	femmodel->parameters->FindParam(&utide, BasalforcingsLaddieVelTideEnum);

	/*Get basal friction coefficient*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieVelFrictionEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		Input* vx_input    = element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
		Input* vy_input    = element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);

		/*friction velocity*/
		IssmDouble* ustar=xNew<IssmDouble>(numvertices);

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			
			/*Prepare value*/
			vx_input->GetInputValue(&vx, gauss);
			vy_input->GetInputValue(&vy, gauss);

			/*Calculate friction velocity*/
			ustar[iv] = sqrt(Cdtop*(vx*vx + vy*vy + utide*utide));
		}

		/*Assign value in: */
		element->AddInput(BasalforcingsLaddieVelFrictionEnum,ustar,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(ustar);
		delete gauss;
	}
}/*}}}*/
void UpdateLaddieHeatAndSaltExchangeCoefficientx(FemModel* femmodel){/*{{{*/
	/*
	 Update gammaT and gammaS (heat exchange coefficient at each vertices

	 See also
	 Eq. 11 and Eq. 12 in Lambert et al. (2023).
	 */

	/*Initialize parameters
	NOTE: These specific parameters are fixed!
	*/
	IssmDouble Pr=13.8; /*Prandtl number*/
	IssmDouble Sc=2432; /*Schmidt number*/
   IssmDouble nu0=1.95e-6; /*Kinematic viscosity = mu / rho*/

	/*Constant value*/
	IssmDouble Cdtop; /*drag coefficient of friction velocity*/

	/*Input value*/
	IssmDouble D;
	IssmDouble ustar;

	femmodel->parameters->FindParam(&Cdtop, BasalforcingsLaddieCdTopEnum);

	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieGammaTEnum,values,P1DGEnum);
			element->AddInput(BasalforcingsLaddieGammaSEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		Input* D_input     = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(D_input);
		Input* ustar_input = element->GetInput(BasalforcingsLaddieVelFrictionEnum); _assert_(ustar_input);

		IssmDouble* gammaT=xNew<IssmDouble>(numvertices);
		IssmDouble* gammaS=xNew<IssmDouble>(numvertices);

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			D_input->GetInputValue(&D,gauss);
			ustar_input->GetInputValue(&ustar,gauss);

			/*Calculate exchange coefficient for temperature*/
			gammaT[iv] = ustar/(2.12*log10(ustar*D/nu0) + 12.5*pow(Pr,2/3) - 8.68);
			gammaS[iv] = ustar/(2.12*log10(ustar*D/nu0) + 12.5*pow(Sc,2/3) - 8.68);
		}

		/*Assign values*/
		element->AddInput(BasalforcingsLaddieGammaTEnum,gammaT,P1DGEnum);
		element->AddInput(BasalforcingsLaddieGammaSEnum,gammaS,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(gammaT);
		xDelete<IssmDouble>(gammaS);
		delete gauss;
	}
}/*}}}*/

void UpdateLaddieMeltratex(FemModel* femmodel){/*{{{*/
	/*
	 Update sub-ice shelf melting from LADDIE depending on melting scheme
	 */

	int ismelt;

	femmodel->parameters->FindParam(&ismelt, BasalforcingsLaddieIsMeltEnum);
	switch(ismelt){
		case 0: /*two equation formulation*/
			_error_("Given md.basalforcings.ismelt=0 is not implemented yet.");
			LaddieMeltrateTwoEquationx(femmodel);
			break;
		case 1: /*three equation formulation*/
			LaddieMeltrateThreeEquationx(femmodel);
			break;
		default:
			_error_("Not implemented yet.");
	}
}/*}}}*/
void LaddieMeltrateTwoEquationx(FemModel *femmodel){/*{{{*/
	/*
	Calculate melt rate with two equation formulation from Holland et al. (2006).
	??
	 */
	_error_("Not implemented yet.");
}/*}}}*/
void LaddieMeltrateThreeEquationx(FemModel* femmodel){ /*{{{*/
	/*
	 Calculate sub-ice shelf melting rate using three equation formulation.

	 See also
	 laddie > src > physics.py > def update_melt(object)
	 */

	/*Initialize input values: */
	int numvertices;
	bool        isgammaTfix;
	IssmDouble  D, vx, vy, T, S, zb;
	IssmDouble  That, Chat, Ctil;
	IssmDouble  AA;
	IssmDouble  b, c;
	IssmDouble  gammaT_const;

	/*Initialize constant values: */
	IssmDouble Utide; /*Tidal velocity [m s-1]*/
   IssmDouble nu0=1.95e-6; /*Kinematic viscosity = mu / rho*/
	IssmDouble Pr=13.8; /*Prandtl number*/
	IssmDouble Sc=2432; /*Schmidt number*/
	IssmDouble Cd_top; /*top drag coefficient*/

	IssmDouble Ci; /*Heat capacity of ice*/
	IssmDouble Cp; /*Heat capacity of seawater*/
	IssmDouble L; /*Latent heat of fusion*/

	IssmDouble l1=-5.73e-2; /*PMP salinity parameter [degC psu-1]*/
	IssmDouble l2=8.32e-2; /*PMP salinity parameter [degC]*/
	IssmDouble l3=7.61e-4; /*PMP salinity parameter [degC m-1]*/

	IssmDouble Ti; /*basal ice temperature on ice shelf [degC]*/

	IssmDouble temp; /*Dummy value*/

	IssmDouble *values; /*value for assigning value*/
	IssmDouble *melt, *Tb;
	IssmDouble  ustar;
	IssmDouble *gammaT, *gammaS;
	Input   *D_input, *vx_input, *vy_input, *T_input, *S_input, *base_input;
	Input   *ustar_input;
	Element *element;
	Gauss   *gauss;

	/*Get parameters*/
	femmodel->parameters->FindParam(&Utide,  BasalforcingsLaddieVelTideEnum);
	femmodel->parameters->FindParam(&Ti,  BasalforcingsLaddieIceTemperatureEnum);
	femmodel->parameters->FindParam(&Cd_top, BasalforcingsLaddieCdTopEnum);
	femmodel->parameters->FindParam(&Ci, MaterialsHeatcapacityEnum);
	femmodel->parameters->FindParam(&Cp, MaterialsMixedLayerCapacityEnum);
	femmodel->parameters->FindParam(&L, MaterialsLatentheatEnum);
	femmodel->parameters->FindParam(&isgammaTfix, BasalforcingsLaddieIsGammaTFixEnum);
	femmodel->parameters->FindParam(&gammaT_const,BasalforcingsLaddieConstGammaTEnum);

	/*Update sub-ice shelf melting rate using three equation formulatoin*/
	for(Object* & object : femmodel->elements->objects){
		element = xDynamicCast<Element*>(object);
		numvertices=element->GetNumberOfVertices();

		if (!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieMeltingRateEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);

			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieTbEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);

			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieGammaTEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);

			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieGammaSEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}
		
		/*Get all input*/
		D_input    = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(D_input);
		vx_input   = element->GetInput(BasalforcingsLaddieVxEnum);    _assert_(vx_input);
		vy_input   = element->GetInput(BasalforcingsLaddieVyEnum);    _assert_(vy_input);
		T_input    = element->GetInput(BasalforcingsLaddieTEnum);     _assert_(T_input);
		S_input    = element->GetInput(BasalforcingsLaddieSEnum);     _assert_(S_input);
		base_input = element->GetInput(BaseEnum);                     _assert_(base_input);
		ustar_input= element->GetInput(BasalforcingsLaddieVelFrictionEnum); _assert_(ustar_input);

		melt  = xNew<IssmDouble>(numvertices);
		Tb    = xNew<IssmDouble>(numvertices); /* freezing temperature at ice shelf base*/
		gammaT= xNew<IssmDouble>(numvertices);
		gammaS= xNew<IssmDouble>(numvertices);

		gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			
			/*Prepare value*/
			D_input->GetInputValue(&D, gauss);
			vx_input->GetInputValue(&vx, gauss);
			vy_input->GetInputValue(&vy, gauss);
			T_input->GetInputValue(&T, gauss);
			S_input->GetInputValue(&S, gauss);
			base_input->GetInputValue(&zb, gauss);
			ustar_input->GetInputValue(&ustar, gauss);

			/*Calculate effective gammaT and gammaS*/
			if(isgammaTfix){
				/*
				Use fixed heat exchange coefficient described in Laddie simulation.
				
				Ratio gammaT/gammaS from Asay-David (ISOMIP) experiment
				*/
				gammaT[iv] = gammaT_const;
				gammaS[iv] = gammaT_const/10.0;
			}else{
				AA = 2.12*log(ustar*D/nu0+1e-12); 
				gammaT[iv] = ustar/(AA + 12.5*pow(Pr,2/3) - 8.68);
				gammaS[iv] = ustar/(AA + 12.5*pow(Sc,2/3) - 8.68);
			}

			/*Solve three equation*/
			That = (l2 + l3*zb);
			Chat = Cp/(L - Ci*Ti);
			Ctil = Ci/Cp;

			b = Chat*gammaT[iv]*(That - T) + gammaS[iv]*(1+Chat*Ctil*(That + l1*S));
			c = Chat*gammaT[iv]*gammaS[iv]*(That-T+l1*S);

			/*Melt rate*/
			melt[iv]=0.5*(-b + sqrt(pow(b,2.0) - 4.0*c));
			//temp = b*b - 4.0*c;
			//if (temp < 0.0){
			//	melt[iv]=0.0;
			//}else{
			//	melt[iv]=0.5*(-b + sqrt(pow(b,2.0) - 4.0*c));
			//}

			/*Temperature at ice shelf base*/
			Tb[iv] = (Chat*gammaT[iv]*T - melt[iv])/(Chat*gammaT[iv] + Chat*Ctil*melt[iv]);
			if (xIsNan<IssmDouble>(Tb[iv])){
				_error_("Error: Tb[" << iv << "]" << " encounters Nan value.\n" << 
						 "See below values\n" <<  
						 "b    = " << b << "\n" << 
						 "c    = " << c << "\n" << 
						 "That = " << That << "\n" << 
						 "Chat = " << Chat << "\n" << 
						 "Ctil = " << Ctil << "\n" <<  
						 "gammaT= " << gammaT[iv] << "\n" << 
						 "gammaS= " << gammaS[iv] << "\n" << 
						 "melt= " << melt[iv] << "\n"
						 "Tb  = " << Tb[iv] << "\n"
						 );
			}
			if (xIsInf<IssmDouble>(Tb[iv])){
				_error_("Error: Tb[" << iv << "]" << " encounters Inf value.\n");
			}
		}

		/*Assigne value*/
		element->AddInput(BasalforcingsLaddieMeltingRateEnum,melt,P1DGEnum);
		element->AddInput(BasalforcingsLaddieTbEnum,Tb,P1DGEnum);
		element->AddInput(BasalforcingsLaddieGammaTEnum,gammaT,P1DGEnum);
		element->AddInput(BasalforcingsLaddieGammaSEnum,gammaS,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(melt);
		xDelete<IssmDouble>(Tb);
		xDelete<IssmDouble>(gammaT);
		xDelete<IssmDouble>(gammaS);
		delete gauss;
	}
}/*}}}*/

void UpdateLaddieEntrainmentRatex(FemModel* femmodel){/*{{{*/
	/*
	 Calculate entrainment rate!
	 */
	int numvertices;
	int isentrainment;

	/*Parameters*/
	IssmDouble  mu; /*Parameter in Gaspar parameterization. Gaspar: 0.5; Gladish: 2.5*/
	IssmDouble  Kparam; /*Kparam: Kochergin entrainment rate*/
	IssmDouble  g; /*gravitational acceleration [m s-1]*/
	IssmDouble  maxdentr; /*maximum detrainment rate [m s-1]*/
	IssmDouble  Dmin; /*minimum plume thickness*/
	/*Parameters for refreezing point*/
	IssmDouble  l1=-5.73e-2;
	IssmDouble  l2=8.32e-2;
	IssmDouble  l3=7.61e-4;
	/*Parameters for thermal expansion and haline contraction coefficients*/
	IssmDouble  alpha=3.733e-5; /*[degC-1]*/
	IssmDouble  beta=7.843e-4; /*[psu-1]*/

	/*Initialize input values*/
	IssmDouble  Kh, Ah;
	IssmDouble  dt;
	IssmDouble  T, S;
	IssmDouble  Tb, Sb;
	IssmDouble  zb;
	IssmDouble  drhoa; /*buoyancy difference between mixed layer and ambient ocean interface.*/
	IssmDouble  drhob; /*buoyancy difference between mixed layer and water in contact with ice*/
	
	IssmDouble  drho;
	IssmDouble  vx, vy, ga, thickness, ustar;
	IssmDouble  dvx[2],dvy[2];
	IssmDouble  dvxdx, dvydy;
	IssmDouble  dthk[2];
	IssmDouble  dthkdx, dthkdy;
	IssmDouble  melt;
	IssmDouble *xyz_list=NULL;
	IssmDouble  Ri, Sc;

	IssmDouble *values;
	IssmDouble *entr;
	IssmDouble *entr_dummy, *entr2_dummy; /*entrainment rate*/
	IssmDouble *dentr; /*dentrainment rate*/
	IssmDouble *convD; /*diverence of (Hv) value*/

	Element *element;
	Gauss   *gauss;

	/*Retrieve all parameters: */
	femmodel->parameters->FindParam(&dt, BasalforcingsLaddieSubTimestepEnum);
	femmodel->parameters->FindParam(&g, ConstantsGEnum);
	femmodel->parameters->FindParam(&isentrainment, BasalforcingsLaddieIsEntrainmentEnum);
	femmodel->parameters->FindParam(&Kparam, BasalforcingsLaddieKparamEnum);
	femmodel->parameters->FindParam(&Ah, BasalforcingsLaddieHorizontalViscosityEnum);
	femmodel->parameters->FindParam(&Kh, BasalforcingsLaddieHorizontalDiffusivityEnum);
	femmodel->parameters->FindParam(&Dmin, BasalforcingsLaddieThicknessMinEnum);
	femmodel->parameters->FindParam(&maxdentr, BasalforcingsLaddieMaxDentrainmentEnum);
	femmodel->parameters->FindParam(&mu, BasalforcingsLaddieMuEnum);

	/*Update entrainment rate*/
	for(Object* & object : femmodel->elements->objects){
		element = xDynamicCast<Element*>(object);
		numvertices=element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieEntrainmentRateEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);

			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieDEntrainmentRateEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);

			values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieConvDEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		/*Retrieve all inputs: */
		Input *ga_input    = element->GetInput(BasalforcingsLaddieAmbientGEnum); _assert_(ga_input);
		Input *drho_input  = element->GetInput(BasalforcingsLaddieDRhoEnum); _assert_(drho_input);

		Input *S_input     = element->GetInput(BasalforcingsLaddieSEnum); _assert_(S_input);

		Input *T_input     = element->GetInput(BasalforcingsLaddieTEnum); _assert_(T_input);
		Input *Tb_input    = element->GetInput(BasalforcingsLaddieTbEnum); _assert_(Tb_input);
		Input *ustar_input = element->GetInput(BasalforcingsLaddieVelFrictionEnum); _assert_(ustar_input);
		Input *zb_input    = element->GetInput(BaseEnum); _assert_(zb_input);

		Input *vx_input=element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
		Input *vy_input=element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);
		Input *thickness_input=element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(thickness_input);
		Input *melt_input=element->GetInput(BasalforcingsLaddieMeltingRateEnum); _assert_(melt_input);

		entr        = xNew<IssmDouble>(numvertices);
		dentr       = xNew<IssmDouble>(numvertices);
		entr_dummy  = xNew<IssmDouble>(numvertices);
		entr2_dummy = xNew<IssmDouble>(numvertices);
		convD       = xNew<IssmDouble>(numvertices);

		gauss=element->NewGauss();
		if(isentrainment==0){/*{{{*/
			for (int iv=0;iv<numvertices;iv++){
				gauss->GaussVertex(iv);

				/*Update entrainment rate with Holland et al. (2006): doi:10.1175/JPO2970.1*/
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				thickness_input->GetInputValue(&thickness,gauss);
				ga_input->GetInputValue(&ga,gauss);
				drho_input->GetInputValue(&drho,gauss);

				if(true){
					Ri = ga*thickness/(pow(vx,2.0) + pow(vy,2.0) + 1e-8);
					Sc = Ri/0.725/(Ri + 0.186 - sqrt(Ri*Ri - 0.316*Ri + 0.0346) + 1e-8);
					
					entr_dummy[iv] = pow(Kparam,2.0)/Sc*sqrt((vx*vx + vy*vy)*(1 + Ri/Sc));
					dentr[iv] = 0.0;
				}else{
					IssmDouble tmp1, tmp2;
					tmp1 = Kparam*Kparam/(Ah*Ah);
					tmp2 = sqrt((vx*vx + vy*vy) - g*drho*Kh/Ah*thickness);

					entr_dummy[iv] = tmp1*tmp2;
					dentr[iv] = 0.0;
				}
			}
		}/*}}}*/
		else if(isentrainment==1){ /*{{{*/
			/*
			Update entrainment rate with Gaspar et al. (1988): doi:10.3189/2012JoG12J003
			
			See also
			--------
			laddie > src > physics.py > "def update_entrainment" in laddie.

			Gaspar et al. (1988): https://doi.org/10.1175/1520-0485(1988)018<0161:MTSCOT>2.0.CO;2.
			Gladish et al. (2012) doi:10.3189/2012JoG12J003
			 */

			for (int iv=0;iv<numvertices;iv++){
				gauss->GaussVertex(iv);

				T_input->GetInputValue(&T,gauss);
				Tb_input->GetInputValue(&Tb,gauss);
				S_input->GetInputValue(&S,gauss);
				zb_input->GetInputValue(&zb,gauss);
				thickness_input->GetInputValue(&thickness,gauss);
				ustar_input->GetInputValue(&ustar,gauss);
				drho_input->GetInputValue(&drhoa,gauss);

				/*Guess salinity under the ice*/
				Sb = (Tb-l2-l3*zb)/l1;
				drhob = (beta*(S-Sb) - alpha*(T-Tb));

				/*Dummy entrainment*/
				entr_dummy[iv] = 2*mu/g*pow(ustar,3.0)/(thickness*drhoa) - drhob/drhoa*melt;
				if (xIsNan<IssmDouble>(entr_dummy[iv]) || xIsInf<IssmDouble>(entr_dummy[iv])){
					_error_("entr_dummy[" << iv << "] got NaN/Inf value!\n" << \
								"mu         = " << mu << "\n" << \
								"g          = " << g << "\n" << \
								"ustar      = " << ustar << "\n" << \
								"thickness  = " << thickness << "\n" << \
								"drhoa      = " << drhoa << "\n" << \
								"drhob      = " << drhob << "\n" << \
								"Tb         = " << Tb << "\n" << \
								"Sb         = " << Sb << "\n" << \
								"T          = " << T << "\n" << \
								"S          = " << S << "\n" << \
								"melt       = " << melt << "\n");
				}
				entr_dummy[iv] = max(entr_dummy[iv], 0.0);
				dentr[iv] = min(maxdentr,max(entr_dummy[iv],0.0));
			}
		}/*}}}*/
		else{/*{{{*/
			_error_("Given md.basalforcings.isentrainment is not avavilable. Only 0 or 1 are available");
		}/*}}}*/

		/*Check Nan value*/
		for(int iv=0;iv<numvertices;iv++){
			if (xIsNan<IssmDouble>(entr_dummy[iv])){
				_error_("entr_dummy[" << iv << "] got NaN value!\n");
			}
			if (xIsInf<IssmDouble>(entr_dummy[iv])){
				_printf0_("entr_dummy value = " << entr_dummy[iv] << "\n");
				_error_("entr_dummy[" << iv << "] got Inf value!\n");
			}
			if (xIsNan<IssmDouble>(dentr[iv])){
				_error_("dentr[" << iv << "] got NaN value!\n");
			}
			if (xIsInf<IssmDouble>(dentr[iv])){
				_error_("dentr[" << iv << "] got Inf value!\n");
			}
		}

		/*Additional entrainment to prevent D < Dmin*/
		element->GetVerticesCoordinates(&xyz_list);
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			thickness_input->GetInputValue(&thickness,gauss);
			melt_input->GetInputValue(&melt,gauss);

			vx_input->GetInputValue(&vx,gauss);
			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vy_input->GetInputValue(&vy,gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			thickness_input->GetInputValue(&thickness,gauss);
			thickness_input->GetInputDerivativeValue(&dthk[0],xyz_list,gauss);

			dvxdx=dvx[0];
			dvydy=dvy[1];
			dthkdx=dthk[0];
			dthkdy=dthk[1];

			/*convD operation in Laddie may be -\nabla \cdot (Hv)??...*/
			convD[iv] = -vx*dthkdx - thickness*dvxdx \
					  -vy*dthkdy - thickness*dvydy;
			entr2_dummy[iv] = max(0.0,
						(Dmin - thickness)/dt - (convD[iv] + melt + entr_dummy[iv] - dentr[iv]));

			/*Get net entrainment rate*/
			entr[iv] = entr_dummy[iv] + entr2_dummy[iv] - dentr[iv];
		}

		/*Assign input*/
		element->AddInput(BasalforcingsLaddieEntrainmentRateEnum,entr,P1DGEnum);
		element->AddInput(BasalforcingsLaddieDEntrainmentRateEnum,dentr,P1DGEnum);
		element->AddInput(BasalforcingsLaddieConvDEnum,convD,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(entr);
		xDelete<IssmDouble>(entr_dummy);
		xDelete<IssmDouble>(entr2_dummy);
		xDelete<IssmDouble>(dentr);
		xDelete<IssmDouble>(convD);
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}
}/*}}}*/
