%Test Name: SquareShelfAdolcTransientControls
md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md.initialization.vx(:)=1.;
md.initialization.vy(:)=1.;
md.geometry.thickness(:)=500-md.mesh.x/10000;
md.geometry.bed =-100-md.mesh.x/1000;
md.geometry.base=-md.geometry.thickness*md.materials.rho_ice/md.materials.rho_water;
md.mask.ocean_levelset=md.geometry.thickness+md.materials.rho_water/md.materials.rho_ice*md.geometry.bed;
pos=find(md.mask.ocean_levelset>=0);
md.geometry.base(pos)=md.geometry.bed(pos);
md.geometry.surface=md.geometry.base+md.geometry.thickness;
md=setflowequation(md,'SSA','all');

%control parameters
md.inversion=adm1qn3inversion(md.inversion);
md.inversion.iscontrol=1;
md.autodiff.isautodiff=1;
md.autodiff.driver='fos_reverse';

md.friction.coefficient(1:md.mesh.numberofvertices,1)=50;
md.friction.coefficient(1:md.mesh.numberofvertices,2)=25;
md.friction.coefficient(md.mesh.numberofvertices+1,1:2)=[0.75,1.25];
min_parameters(1:md.mesh.numberofvertices,1:2)=1;
min_parameters(md.mesh.numberofvertices+1,1:2)=[0.75,1.25];
max_parameters(1:md.mesh.numberofvertices,1:2)=500;
max_parameters(md.mesh.numberofvertices+1,1:2)=[0.75,1.25];
md.autodiff.independents = {independent('name','FrictionCoefficient',...
	'control_size',2,...
	'type','vertex',...
	'min_parameters',min_parameters,...
	'max_parameters',max_parameters,...
	'control_scaling_factor',1)...
	};

md.outputdefinition.definitions{1}=cfsurfacesquare('name','VyMisfit1',...
	'definitionstring','Outputdefinition1',...
	'model_string','Vy',...
	'observation_string','VyObs',...
	'observation',md.initialization.vy/md.constants.yts,...
	'weights',ones(md.mesh.numberofvertices,1),...
	'weights_string','WeightsSurfaceObservation',...
	'datatime',0.75);

md.timestepping.interp_forcing=1;
md.timestepping.time_step=0.5;
md.timestepping.final_time=1.5;

md.transient.ismasstransport=1;
md.transient.isstressbalance=1;
md.transient.isgroundingline=1;
md.transient.ismovingfront=0;
md.transient.isthermal=0;

pos=find(md.mask.ocean_levelset<0);
md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=25*zeros(md.mesh.numberofvertices,1);

md.autodiff.dependents{1} = dependent('name','Outputdefinition1','type','scalar','fos_reverse_index',1);
md.inversion.maxiter = 2;
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'transient');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','FrictionCoefficient','Pressure1','Vel1','Vx1','Vy1','Pressure2','Vel2','Vx2','Vy2'};
field_tolerances={1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11};
field_values={...
	(md.results.TransientSolution(1).Gradient1),...
	md.results.TransientSolution(1).J,...
	(md.results.TransientSolution(1).FrictionCoefficient),...
	(md.results.TransientSolution(1).Pressure),...
	(md.results.TransientSolution(1).Vel),...
	(md.results.TransientSolution(1).Vx),...
	(md.results.TransientSolution(1).Vy),...
	(md.results.TransientSolution(7).Pressure),...
	(md.results.TransientSolution(7).Vel),...
	(md.results.TransientSolution(7).Vx),...
	(md.results.TransientSolution(7).Vy)
};
