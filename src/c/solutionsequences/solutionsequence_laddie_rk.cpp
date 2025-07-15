/*!\file: solutionsequence_fct.cpp
 * \brief: numerical core of flux corrected transport solution
 */

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../analyses/analyses.h"
#define USEPENALTYMETHOD false

void solutionsequence_laddie_rk(FemModel* femmodel){/*{{{*/
	/*NOTE: Runge-Kutta 2 scheme*/

	/*intermediary: */
	IssmDouble           theta,deltat,dmax;
	int                  configuration_type,analysis_type;
	Matrix<IssmDouble>*  K  = NULL;
	Matrix<IssmDouble>*  Mc = NULL;
	Vector<IssmDouble>*  p  = NULL;
	Vector<IssmDouble>*  ug = NULL;
	Vector<IssmDouble>*  uf = NULL;
	BasalforcingsLaddieMomentumAnalysis* manalysis = NULL;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&deltat,BasalforcingsLaddieSubTimestepEnum);
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->UpdateConstraintsx();

	#ifdef _HAVE_PETSC_
	#ifdef _HAVE_CODIPACK_
		_error_("No CoDiPack handling for PETSc and fct");
	#else
	/*Convert matrices to PETSc matrices*/
	Mat D_petsc  = NULL;
	Vec u        = NULL;
	Vec udot     = NULL;
	Mat K_petsc  = NULL; // K->pmatrix->matrix;
	Mat Mc_petsc = NULL; // Mc->pmatrix->matrix;
	Vec p_petsc	 = NULL; // p->pvector->vector;
	//NOTE: Dummy values.
	//Mat K_petsc  = K->pmatrix->matrix;
	//Mat Mc_petsc = Mc->pmatrix->matrix;
	//Vec p_petsc	 = p->pvector->vector;

	/*un0 and un1 represent u^{n} and u^{n+1}*/
	Vector<IssmDouble>* un0  = NULL;
	Vector<IssmDouble>* un1  = NULL;
	Vec u0       = NULL;
	Vec u1       = NULL;
	Vec u2       = NULL;

	Vec Ku       = NULL; /*K u ~ multiply stiffness matrix and vector*/

	Vec r0       = NULL;
	Vec r1       = NULL;

	Vec r0_dot   = NULL;
	Vec r1_dot   = NULL;

	/*Get previous solution u^n*/
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&un0, ug, femmodel->nodes,femmodel->parameters);
	delete ug;

	/*Now, solve RK2...*/

	/*Step1:
	  Solve
	  u0 = un0
	  u1 = un0 + dt * 1 * M-1 * r(u0)
	 */

	/*Recovert parameters: */
	femmodel->UpdateConstraintsx();

	/*Create analysis*/
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	switch(analysis_type){
		case BasalforcingsLaddieMomentumAnalysisEnum:
			manalysis = new BasalforcingsLaddieMomentumAnalysis();
			manalysis->RkMassMatrix(&Mc,femmodel);
			manalysis->RkKMatrix(&K,femmodel);
			manalysis->RkPVector(&p,femmodel);
			break;
		default: _error_("analysis type " << EnumToStringx(analysis_type) << " not supported for FCT\n");
	}
	delete manalysis;

	/*Copy IssmDouble to PetscVec/Mat*/
	_printf0_("   Copy variable IssmDouble to Vec\n");
	u0      =un0->pvector->vector;
	Mc_petsc=Mc->pmatrix->matrix;
	K_petsc =K->pmatrix->matrix;
	r0      =p->pvector->vector;

	/*Go solve solution*/
	_printf0_("RK2 step1 \n");
	MatMult(K_petsc, u0, Ku);
	VecAXPY(r0,1.0,Ku); /*y = ax + y*/

	femmodel->profiler->Start(SOLVER);
	SolverxPetsc(&r0_dot,Mc->pmatrix->matrix,r0,NULL,NULL, femmodel->parameters);
	femmodel->profiler->Stop(SOLVER);

	/*u^1 = u^0 + dt * 1 * M^{-1} * r^{0}*/
	VecWAXPY(u1, 1.0*deltat, r0_dot, u0); /* VecWAXPY: w = alpha * x + y */
	MatFree(&Mc_petsc);
	MatFree(&K_petsc);
	VecFree(&p_petsc);
	VecFree(&Ku);

	/*Update element input*/
	uf=new Vector<IssmDouble>(u1);
	InputUpdateFromSolutionx(femmodel, uf);
	delete uf;

	/*Step2: Final step for RK2
		Solve
		u^{n+1} = u^{n} + dt * (0.5 * M^{-1} r(u0)  + 0.5 * M^{-1} r(u1))
	 
	 */
	/*Recovert parameters: */
	femmodel->UpdateConstraintsx();

	switch(analysis_type){
		case BasalforcingsLaddieMomentumAnalysisEnum:
			manalysis = new BasalforcingsLaddieMomentumAnalysis();
			manalysis->RkMassMatrix(&Mc,femmodel);
			manalysis->RkKMatrix(&K,femmodel);
			manalysis->RkPVector(&p,femmodel);
			break;
		default: _error_("analysis type " << EnumToStringx(analysis_type) << " not supported for FCT\n");
	}
	delete manalysis;

	K_petsc  = K->pmatrix->matrix;
	Mc_petsc = Mc->pmatrix->matrix;
	r1       = p->pvector->vector;

	/*Go solve solution*/
	femmodel->profiler->Start(SOLVER);
	_printf0_("RK2 step2 \n");
	SolverxPetsc(&r1_dot,Mc->pmatrix->matrix,r1,NULL,NULL, femmodel->parameters);
	femmodel->profiler->Stop(SOLVER);

	/* u(n+1)  = dt*0.5*r0_dot + u(n) */
	VecWAXPY(u2, 0.5*deltat, r0_dot, u0);
	/* u(n+1) += dt*0.5*r1_dot*/
	VecAXPY(u2, 0.5*deltat, r1_dot);

	/*Serialize u and udot*/
	//IssmDouble* udot_serial = NULL;
	//IssmDouble* u_serial    = NULL;
	//IssmDouble* ml_serial   = NULL;
	//VecToMPISerial(&udot_serial,udot    ,IssmComm::GetComm());
	//VecToMPISerial(&u_serial   ,u       ,IssmComm::GetComm());
	//VecToMPISerial(&ml_serial  ,Ml_petsc,IssmComm::GetComm());

	/*Clean up*/
	//MatFree(&D_petsc);
	//delete Mc;
	//xDelete<IssmDouble>(udot_serial);
	//xDelete<IssmDouble>(u_serial);
	//xDelete<IssmDouble>(ml_serial);

	///*Compute solution u^n+1 = u_L + deltat Ml^-1 fbar*/
	//UpdateSolution(u,udot,Ml_petsc,Fbar,deltat);
	//uf =new Vector<IssmDouble>(u);
	//VecFree(&u);
	//VecFree(&Fbar);
	//VecFree(&udot);
	//VecFree(&Ml_petsc);
	
	uf = new Vector<IssmDouble>(u2);
	delete un0;
	delete un1;
	VecFree(&u0);
	VecFree(&u1);
	VecFree(&u2);
	VecFree(&r0);
	VecFree(&r1);
	VecFree(&r0_dot);
	VecFree(&r1_dot);

	/*Update Element inputs*/
	InputUpdateFromSolutionx(femmodel,uf);
	delete uf;

	#endif
	#else
	_error_("PETSc needs to be installed");
	#endif
}/*}}}*/
