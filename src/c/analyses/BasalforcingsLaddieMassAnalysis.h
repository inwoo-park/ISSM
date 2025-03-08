/*! \file BasalforcingsPlumeAnalysis.h 
 *  \brief: header file for generic external result object
 */

 #ifndef _BasalforcingsPlumeAnalysis_
 #define _BasalforcingsPlumeAnalysis_
 
 /*Headers*/
 #include "./Analysis.h"
 
 class BasalforcingsPlumeAnalysis: public Analysis{
 
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

			/*Mass balance special*/
			ElementMatrix* CreateKMatrixMassbalance(Element* element);
			ElementMatrix* CreateKMatrixMassbalanceCG(Element* element);
			ElementVector* CreatePVectorMassbalance(Element* element);
			ElementVector* CreatePVectorMassbalanceCG(Element* element);

         void           GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element);
         void           GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index);
         void           InputUpdateFromSolution(IssmDouble* solution,Element* element);
         void           UpdateConstraints(FemModel* femmodel);
 };
 #endif
 
/*External functions*/
IssmDouble GetBasalforcingsFrictionVelocity(IssmDouble Cd_top, IssmDouble vx, IssmDouble vy, IssmDouble Utide);
void GetHeatExchangeCoefficient(IssmDouble *pgammaT, IssmDouble *pgammaS, IssmDouble Ustar, IssmDouble D);
IssmDouble GetEffectiveGravitationDensity(IssmDouble g, IssmDouble Ta, IssmDouble Sa, IssmDouble T, IssmDouble S);

IssmDouble GetEntrainmentRate(IssmDouble isentrainment, IssmDouble ga, IssmDouble gb, IssmDouble depth, IssmDouble vx, IssmDouble vy, IssmDouble meltrate);
