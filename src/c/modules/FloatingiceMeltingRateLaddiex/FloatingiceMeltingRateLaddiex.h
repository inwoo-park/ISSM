/*!\file:  FloatingiceMeltingRateLaddiex.h
 * \brief header file for external functions for LADDIE model
 */ 

#ifndef _FloatingiceMeltingRateLaddiex_H
#define _FloatingiceMeltingRateLaddiex_H

#include "../../classes/classes.h"
#include "../FloatingiceMeltingRatex/FloatingiceMeltingRatex.h"

IssmDouble GetDensityDifferencex(IssmDouble rho0, IssmDouble T, IssmDouble S, IssmDouble Ta, IssmDouble Sa);
IssmDouble GetEffectiveGravitationAccelerationx(IssmDouble g, IssmDouble rho0, IssmDouble T, IssmDouble S, IssmDouble Ta, IssmDouble Sa);
void FloatingiceMeltingRateLaddiex(FemModel* femmodel);

void UpdateLaddieDThicknessDtx(FemModel* femmodel);
void UpdateLaddieAmbientFieldx(FemModel* femmodel);
void UpdateLaddieDensityAndEffectiveGravityx(FemModel* femmodel);

void UpdateLaddieFrictionVelocityx(FemModel* femmodel);
void UpdateLaddieHeatAndSaltExchangeCoefficientx(FemModel* femmodel);

void UpdateLaddieMeltratex(FemModel* femmodel);
void LaddieMeltrateTwoEquationx(FemModel *femmodel);
void LaddieMeltrateThreeEquationx(FemModel* femmodel);

void UpdateLaddieEntrainmentRatex(FemModel* femmodel);
IssmDouble GetEntrainmentRateHollandx(IssmDouble Kparam, IssmDouble ga, IssmDouble thickness, IssmDouble vx, IssmDouble vy);
IssmDouble GetEntrainmentRateGasparx(IssmDouble g, IssmDouble rho0, IssmDouble D, IssmDouble T, IssmDouble S, IssmDouble Ta, IssmDouble Sb, IssmDouble Ustar, IssmDouble melt_rate);

#endif
