/*!\file FloatingiceMeltingRateLaddiex
 * \brief: calculates Floating ice melting rate
 */

#include "./FloatingiceMeltingRateLaddiex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "./../../classes/Inputs/DatasetInput.h"
#include "../InputDuplicatex/InputDuplicatex.h"

IssmDouble GetLaddieFrictionVelocityx(IssmDouble Cd_top, IssmDouble vx, IssmDouble vy, IssmDouble Utide){/*{{{*/
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
	 Eq 13. in Lambert et al. (2023TC)
	 */

	/*Hard coding for specific parameters*/
	//IssmDouble Cd_top=1.1e-3; // top drag coefficient 
	//IssmDouble Utide=0.01; // tide velocity [m s-1]
	
	/*Input values*/
	//IssmDouble vx, vy;

	/*Output values*/
	IssmDouble Ustar; // Friction velocity

	/*Calculate fricition velocity*/
	Ustar = Cd_top * (pow(vx,2.0) + pow(vy,2.0) + pow(Utide,2));
   Ustar = pow(Ustar,0.5);

	return Ustar;
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
void FloatingiceMeltingRateLaddiex(FemModel* femmodel){/*{{{*/
	int ismelt;
	femmodel->parameters->FindParam(&ismelt, BasalforcingsLaddieVelTideEnum);

	switch(ismelt){
		case 0: /*Two equation formulation*/
			_error_("Two equation formulation (ismelt = 0) is not implemented yet.");
			break;
		case 1: /*Three equation formulation*/
			LaddieMeltrateThreeEquationx(femmodel);
			break;
	}
}/*}}}*/

void UpdateLaddieAmbientFieldx(FemModel* femmodel){/*{{{*/
	/*
	Update ambient ocean temperature and salinity from forcing fields.

	Make 3D ambient ocean temperature/salinity to 2D ambient ocean temperature/salinity considering ice base elevation.
	 */
	
	IssmDouble* forcing_depth=NULL;
	int num_depths;

	femmodel->parameters->FindParam(&forcing_depth,&num_depths,BasalforcingsLaddieForcingDepthEnum); _assert_(depths);

	/*Get ambient ocean temperature and salinity at each ice shelf point - linearly interpolate in depth and time*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieAmbientTemperatureEnum,values,P1DGEnum);
			element->AddInput(BasalforcingsLaddieAmbientSalinityEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}
		
		/*Get ambient ocean temeprature and salinity on all vertices*/
		IssmDouble* Ttmp           = xNew<IssmDouble>(numvertices);
		IssmDouble* Stmp           = xNew<IssmDouble>(numvertices);
		IssmDouble* depth_vertices = xNew<IssmDouble>(numvertices);
		/*
		 * tf: temperature forcing
		 * sf: salinity forcing
		 */
		DatasetInput* tf_input = element->GetDatasetInput(BasalforcingsLaddieForcingTemperatureEnum); _assert_(tf_input);
		DatasetInput* sf_input = element->GetDatasetInput(BasalforcingsLaddieForcingTemperatureEnum); _assert_(tf_input);

		element->GetInputListOnVertices(&depth_vertices[0],BaseEnum);

		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			IssmDouble depth=-depth_vertices[iv];
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

			element->AddInput(BasalforcingsLaddieAmbientTemperatureEnum,Ttmp,P1DGEnum);
			element->AddInput(BasalforcingsLaddieAmbientSalinityEnum,Stmp,P1DGEnum);
			xDelete<IssmDouble>(Ttmp);
			xDelete<IssmDouble>(Stmp);
			xDelete<IssmDouble>(depth_vertices);
			delete gauss;
		}
	}
}/*}}}*/
void UpdateLaddieDensityAndEffectiveGravityx(FemModel* femmodel){/*{{{*/
	/*
	Update ocean density and effective gravitational acceleration due to buoyancy.
	
	See also
	Eq. 6 and Eq. 7 in Lambert et al. (2023): https://tc.copernicus.org/articles/17/3203/2023/
	 */

	IssmDouble rho0; /*default ocean density [kg m-3]*/
	IssmDouble g; /*gravitaional acceleration [m s-1]*/
	IssmDouble T,S; /*plume temperature and salinity*/
	IssmDouble Ta, Sa; /*ambient ocean forcing for temperature [degC] and salinity [psu]*/
	IssmDouble alpha=3.733e-5; /*thermal expansion coefficient [degC -1]*/
	IssmDouble beta=7.843e-4; /*haline contraction coefficient [psu-1]*/

	femmodel->parameters->FindParam(&rho0,MaterialsRhoSeawaterEnum);
	femmodel->parameters->FindParam(&g,ConstantsGEnum);

	/*Get basal friction coefficient*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieDRhoEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		/*Get inputs*/
		/*Plume temperature and salinity*/
		Input* T_input = element->GetInput(BasalforcingsLaddieTEnum); _assert_(T_input);
		Input* S_input = element->GetInput(BasalforcingsLaddieSEnum); _assert_(S_input);
		/*Ambient ocean temperature and salinity*/
		Input* Ta_input = element->GetInput(BasalforcingsLaddieAmbientTemperatureEnum); _assert_(Ta_input);
		Input* Sa_input = element->GetInput(BasalforcingsLaddieAmbientTemperatureEnum); _assert_(Sa_input);

		/*Initialize variable*/
		IssmDouble* ga=xNew<IssmDouble>(numvertices);
		IssmDouble* drho=xNew<IssmDouble>(numvertices);

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			
			/*Prepare value*/
			T_input->GetInputValue(&T, gauss);
			S_input->GetInputValue(&S, gauss);
			Ta_input->GetInputValue(&Ta, gauss);
			Sa_input->GetInputValue(&Sa, gauss);

			/*Calculate density difference*/
			drho[iv] = rho0*(-alpha*(Ta-T) + beta*(Sa-S));
			ga[iv] = g*drho[iv]/rho0;
		}

		/*Assign value in: */
		element->AddInput(BasalforcingsLaddieDRhoEnum,drho,P1DGEnum);
		element->AddInput(BasalforcingsLaddieAmbientGEnum,drho,P1DGEnum);
		
		/*Clear memory*/
		xDelete<IssmDouble>(ga);
		xDelete<IssmDouble>(drho);
	}
}/*}}}*/

void UpdateLaddieFrictionVelocityx(FemModel* femmodel){/*{{{*/
	/*
	 Update friction velocity, Ustar, at each vertices.

	 See also
	 Eq. 13 in Lambert et al. (2023): https://tc.copernicus.org/articles/17/3203/2023/
	 */
	IssmDouble Cdtop; /*Friction coefficient [m s-1]*/
	IssmDouble vx, vy; /*plume velocity vector [m s-1]*/
	IssmDouble utide; /*tide velocity*/

	femmodel->parameters->FindParam(&Cdtop, BasalforcingsLaddieCdTopEnum);

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
		Input* utide_input = element->GetInput(BasalforcingsLaddieVelTideEnum); _assert_(utide_input);

		/*friction velocity*/
		IssmDouble* ustar=xNew<IssmDouble>(numvertices);

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			
			/*Prepare value*/
			vx_input->GetInputValue(&vx, gauss);
			vy_input->GetInputValue(&vy, gauss);
			utide_input->GetInputValue(&utide, gauss);

			/*Calculate friction velocity*/
			ustar[iv] = GetLaddieFrictionVelocityx(Cdtop, vx, vy, utide);
		}

		/*Assign value in: */
		element->AddInput(BasalforcingsLaddieVelFrictionEnum,ustar,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(ustar);
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
		Input* ustar_input = element->GetInput(BasalforcingsLaddieVelFrictionEnum); _assert_(ustart_input);

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
		case 1: /*three equation formulation*/
			LaddieMeltrateThreeEquationx(femmodel);
		default:
			_error_("Not implemented yet.");
	}
}/*}}}*/
void LaddieMeltrateTwoEquationx(FemModel *femmodel){/*{{{*/
	/*
	Calculate melt rate with two equation formulation.

	See also
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

	IssmDouble Utide; /*Tidal velocity [m s-1]*/
	IssmDouble Ustar; /*Friction velocity [m s-1]*/
   IssmDouble nu0=1.95e-6; /*Kinematic viscosity = mu / rho*/
	IssmDouble Pr=13.8; /*Prandtl number*/
	IssmDouble Sc=2432; /*Schmidt number*/
	IssmDouble Cd_top; /*top drag coefficient*/

	IssmDouble That, Chat, Ctil;
	IssmDouble Ci; /*Heat capacity of ice*/
	IssmDouble Cp; /*Heat capacity of seawater*/
	IssmDouble L; /*Latent heat of fusion*/

	IssmDouble l1=-5.73e-2; /*PMP salinity parameter [degC psu-1]*/
	IssmDouble l2=8.32e-2; /*PMP salinity parameter [degC]*/
	IssmDouble l3=7.61e-4; /*PMP salinity parameter [degC m-1]*/

	IssmDouble Ti=-20; /*basal ice temperature on ice shelf [degC]*/

	/*Get parameters*/
	femmodel->parameters->FindParam(&Utide,  BasalforcingsLaddieVelTideEnum);
	femmodel->parameters->FindParam(&Cd_top, BasalforcingsLaddieCdTopEnum);
	femmodel->parameters->FindParam(&Ci, MaterialsHeatcapacityEnum);
	femmodel->parameters->FindParam(&Cp, MaterialsMixedLayerCapacityEnum);
	femmodel->parameters->FindParam(&L, MaterialsLatentheatEnum);

	/*Update sub-ice shelf melting rate using three equation formulatoin*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int numvertices=element->GetNumberOfVertices();

		if (!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* value = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsFloatingiceMeltingRateEnum,value,P1DGEnum);
			xDelete<IssmDouble>(value);
			continue;
		}
		
		/*Get all input*/
		Input* D_input    = element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(D_input);
		Input* vx_input   = element->GetInput(BasalforcingsLaddieVxEnum);    _assert_(vx_input);
		Input* vy_input   = element->GetInput(BasalforcingsLaddieVyEnum);    _assert_(vy_input);
		Input* T_input    = element->GetInput(BasalforcingsLaddieTEnum);     _assert_(T_input);
		Input* S_input    = element->GetInput(BasalforcingsLaddieSEnum);     _assert_(S_input);
		Input* base_input = element->GetInput(BaseEnum);                     _assert_(base_input);

		IssmDouble* melt = xNew<IssmDouble>(numvertices);
		IssmDouble  D, vx, vy, T, S, zb;
		IssmDouble  That, Chat, Ctil;
		IssmDouble  AA;
		IssmDouble  b, c;
		IssmDouble  gammaT, gammaS; /*Heat exchange coefficient for three equation*/

		Gauss* gauss=element->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			
			/*Prepare value*/
			D_input->GetInputValue(&D, gauss);
			vx_input->GetInputValue(&vx, gauss);
			vy_input->GetInputValue(&vy, gauss);
			T_input->GetInputValue(&T, gauss);
			S_input->GetInputValue(&S, gauss);
			base_input->GetInputValue(&zb, gauss);

			/*Calculate friction velocity*/
			Ustar = GetLaddieFrictionVelocityx(Cd_top, vx, vy, Utide);

			/*Calculate effective gammaT and gammaS*/
			AA = 2.12*log(Ustar*D/nu0+1e-12); 
			gammaT = Ustar/(AA + 12.5*pow(Pr,2/3) - 8.68);
			gammaS = Ustar/(AA + 12.5*pow(Sc,2/3) - 8.68);

			/*Solve three equation*/
			That = (l2 + l3*zb);
			Chat = Cp/(L - Ci*Ti);
			Ctil = Ci/Cp;

			b = Chat*gammaT*(That - T) + gammaS*(1+Chat*Ctil*(That + l1*S));
			c = Chat*gammaT*gammaS*(That-T+l1*S);

			/*Melt rate*/
			melt[iv]=0.5*(-b + sqrt(pow(b,2.0) - 4*c));
		}

		element->AddInput(BasalforcingsFloatingiceMeltingRateEnum,melt,P1DGEnum);
		xDelete<IssmDouble>(melt);
		delete gauss;
	}
}/*}}}*/

void UpdateLaddieEntrainmentRatex(FemModel* femmodel){/*{{{*/
	/*
	 Calculate entrainment rate!
	 */
	IssmDouble Kparam; /*Kparam: Kochergin entrainment rate*/
	int isentrainment;

	femmodel->parameters->FindParam(&isentrainment, BasalforcingsLaddieIsEntrainmentEnum);

	/*Update sub-ice shelf melting rate using three equation formulatoin*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int numvertices=element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsLaddieEntrainmentRateEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		/*Initialize array*/
		IssmDouble* entr=xNew<IssmDouble>(numvertices);

		switch(isentrainment){
			case 0:

				/*Initialize inputs: */
				Input* vx_input=element->GetInput(BasalforcingsLaddieVxEnum); _assert_(vx_input);
				Input* vy_input=element->GetInput(BasalforcingsLaddieVyEnum); _assert_(vy_input);
				Input* ga_input=element->GetInput(BasalforcingsLaddieAmbientGEnum); _assert_(ga_input);
				Input* D_input=element->GetInput(BasalforcingsLaddieThicknessEnum); _assert_(D_input);

				/*Initialize input values*/
				IssmDouble vx, vy, ga, D;

				Gauss* gauss=element->NewGauss();
				for (int iv=0;iv<numvertices;iv++){
					gauss->GaussVertex(iv);

					vx_input->GetInputValue(&vx,gauss);
					vy_input->GetInputValue(&vy,gauss);
					ga_input->GetInputValue(&ga,gauss);
					D_input->GetInputValue(&D,gauss);
					
					entr[iv]=GetEntrainmentRateHollandx(Kparm, ga, D, vx, vy);
				}
				break;
			case 1:
				_error_("Given md.basalforcings.isentrainement=1 is not implemented yet!");
				break;
			default: _error_("Given md.basalforcings.isentrainment is not avavilable. Only 0 or 1 are available");
		}

		/*Assign input*/
		element->AddInput(BasalforcingsLaddieEntrainmentRateEnum,entr,P1DGEnum);

		/*Clear memory*/
		xDelete<IssmDouble>(entr);
	}
}/*}}}*/
IssmDouble GetEntrainmentRateHollandx(IssmDouble Kparam, IssmDouble ga, IssmDouble thickness, IssmDouble vx, IssmDouble vy){/*{{{*/
	/*
	 Calculate the entrainment rate using Eq. 2 in Holland et al. (2006).
	
	 Inputs
	 * Kparam: Kochergin entrainment rate
	 * ga: effective gravitational acceleration [m s-2]
	 * thickness: plume thickness [m]
	 * u,v: depth averaged plume velocity [m s-1]

	 See also
	 Payne, A. J., Holland, P. R., Shepherd, A. P., Rutt, I. C., Jenkins, A., and Joughin, I.: Numerical modeling of ocean-ice interactions under Pine Island Bay’s ice shelf, J. Geophys. Res., 112, C10019, https://doi.org/10.1029/2006JC003733, 2007.
	 Holland, P. R. and Feltham, D. L.: The Effects of Rotation and Ice Shelf Topography on Frazil-Laden Ice Shelf Water Plumes, Journal of Physical Oceanography, 36, 2312–2327, https://doi.org/10.1175/JPO2970.1, 2006.
	 */
	IssmDouble Ri; /*Richardson number*/
	IssmDouble Sc; /*Schmidt number*/
	IssmDouble entrainment; /*Entrainment rate*/

	/*Calculate Richardson and Schmidt number*/
	Ri = ga*thickness/(pow(vx,2.0) + pow(vy,2.0));
	Sc = Ri/(0.725*(Ri + 0.186 - pow(pow(Ri,2.0) - 0.316*Ri + 0.0346, 0.5)));

	entrainment = pow(Kparam,2.0)/Sc*pow((pow(vx,2.0) + pow(vy,2.0))*(1 + Ri/Sc), 0.5);

	return entrainment;
}/*}}}*/
IssmDouble GetEntrainmentRateGasparx(IssmDouble g, IssmDouble rho0, IssmDouble D, IssmDouble T, IssmDouble S, IssmDouble Ta, IssmDouble Sa, IssmDouble Ustar, IssmDouble melt_rate){/*{{{*/
	/*
	 Calculate the entrainment rate using Eq. 2 in Holland et al. (2006).
	
	 Inputs
	 * g: gravitational acceleration [m s-1]
	 * T, S: ocean temperature and salinity in mixed-layer.
	 * Ta, Sa: ambient ocean temperature and salinity in contacting ocean.
	 * D: plume thickness [m]
	 * Ustar: frictional stress at the ice base [m s-1]
	 * zb: basal elevation [m]
	 * melt_rate: sub-ice shelf melting rate [m s-1]

	 See also
	 Gladish, C. V., Holland, D. M., Holland, P. R., and Price, S. F.: Ice-shelf basal channels in a coupled ice/ocean model, Journal of Glaciology, 58, 1227–1244, https://doi.org/10.3189/2012JoG12J003, 2012.
	 */

	/*Hard coding for specific parameters*/
	IssmDouble mu=0.5; /*Detrainment parameter*/
	IssmDouble ga;
	IssmDouble entrainment; /*Entrainment rate*/

	/*Get effective gravitational acceleration due to density difference between mixed-layer and ambient ocean below the mixed-layer.*/
	ga = GetEffectiveGravitationAccelerationx(g, rho0, T, S, Ta, Sa);

	entrainment = 2*mu*pow(Ustar,3.0)/(ga*D) - melt_rate;

	return entrainment;
}/*}}}*/
