/*! \file BasalforcingsLaddieMomentumAnalysis.h 
 *  \brief: header file for generic external result object
 */

#ifndef _BasalforcingsLaddieMomentumAnalysis_
#define _BasalforcingsLaddieMomentumAnalysis_

/*Headers*/
#include "./Analysis.h"

class BasalforcingsLaddieMomentumAnalysis: public Analysis{

	public:
		/*Model processing*/
      void CreateConstraints(Constraints* constraints,IoModel* iomodel);
      void CreateLoads(Loads* loads, IoModel* iomodel);
      void CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr=false);
      int  DofsPerNode(int** doflist,int domaintype,int approximation);
      void UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type);
      void UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum);
 
      /*Finite element Analysis*/
      void           Core(FemModel* femmodel);
      void           PreCore(FemModel* femmodel);
      ElementVector* CreateDVector(Element* element);
      ElementMatrix* CreateJacobianMatrix(Element* element);
      ElementMatrix* CreateKMatrix(Element* element);
      ElementVector* CreatePVector(Element* element);
      void           GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element);
      void           GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index);
      void           InputUpdateFromSolution(IssmDouble* solution,Element* element);
      void           UpdateConstraints(FemModel* femmodel);

		/*Momentum: CG method*/
		ElementMatrix* CreateKMatrixCG(Element* element);
		ElementVector* CreatePVectorCG(Element* element);
		ElementVector* CreatePVectorCGLinearize(Element* element);

		/*Flux correction transport (FCT)*/
		ElementMatrix* CreateFctKMatrix(Element* element);
		ElementMatrix* CreateMassMatrix(Element* element);
		ElementVector* CreateFctPVector(Element* element);


		void           FctKMatrix(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,FemModel* femmodel);
		void           FctPVector(Vector<IssmDouble>** ppf,FemModel* femmodel);
		void           MassMatrix(Matrix<IssmDouble>** pMff,FemModel* femmodel);

		/*Runge-Kutta scheme*/
		ElementMatrix* CreateKMatrixRk(Element* element);
		ElementVector* CreatePVectorRk(Element* element);
		void           RkKMatrix(Matrix<IssmDouble>** pKff,FemModel* femmodel);
		void           RkMassMatrix(Matrix<IssmDouble>** pMff,FemModel* femmodel);
		void           RkPVector(Vector<IssmDouble>** ppf,FemModel* femmodel);
};
#endif
