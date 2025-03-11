%BASAL FORCINGS class definition
%
%   Usage:
%      basalforcings=basalforcingsladdie();

classdef basalforcingsladdie
	properties (SetAccess=public) % {{{
		%Default parameters
		groundedice_melting_rate  = NaN;
		geothermalflux            = NaN;
		D = NaN;
		vx= NaN;
		vy= NaN;
		temperature = NaN;
		salinity = NaN;

		%Ocean forcing values
		forcing_depth = NaN;
		forcing_temeprature = NaN;
		forcing_salinity = NaN;

		%Dirichlet boundary
		spcD  = NaN;
		spcvx = NaN;
		spcvy = NaN;
		%spctemperature = NaN;
		%spcsalinity    = NaN;

		%Neumann boundary

		%Parameters
		Kh = 0;
		Ah = 0;
		Dmin = 0;
		Utide= 0;
		Cd= 0;
		Cd_top=0;
		isentrainment = 0;
		ismelt=0;

		%Timestepping
		subtimestep=0;
		diagnostic_frequency=0;

		%Numerical Method
		stabilization=0;
	end % }}}
	methods
		function self = basalforcingsladdie(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self =structtoobj(basalforcings(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   basal forcings plume model parameters:'));

			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]');
			fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]');

			fielddisplay(self,'D','initial plume thickness [m].');
			fielddisplay(self,'vx','initial depth averaged plume velocity in x-direction [m s-1].');
			fielddisplay(self,'vy','initial depth averaged plume velocity in y-direction [m s-1].');
			fielddisplay(self,'temperature','initial plume temperature [K].');
			fielddisplay(self,'salinity','initial plume salinity [psu].');

			%parameters
			fielddisplay(self,'Kh','horizontal diffusivity [unit: m2 s-1]');
			fielddisplay(self,'Ah','horizontal viscosity [unit: m2 s-1]');
			fielddisplay(self,'Dmin','minimum plume thickness [unit: m]');
			fielddisplay(self,'Utide','tidal velocity [unit: m s-1]');

			fielddisplay(self,'f_cori','Coriolis frequency [unit: s-1]');
			fielddisplay(self,'Cd','momentum drag coefficient [unit: -]')
			fielddisplay(self,'Cd_top','top drag coefficient [unit: -]')
			fielddisplay(self,'isentrainment','select calculating entrainment value (dot{e}). 0: Holland et al. (2006), 1: Gaspar et al. (1988). (defalt: 0)')
			fielddisplay(self,'ismelt','select calculating sub-ice shelf melting value (M_b). 0: two-equation formulation (McPhee et al., 2008), 1: three-equation formulation (Jenkins et al., 2010). (defalt: 0)')

			fielddisplay(self,'subtimestep','Sub timestepping for 2D plume model [unit: s]');
			fielddisplay(self,'diagnostic_frequency','Store the estimated sub-ice shelf melting rate every given diagnostic_frequency [unit: -]');
		end % }}}
		function self = extrude(self,md) % {{{
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1); 
			self.floatingice_melting_rate=project3d(md,'vector',self.floatingice_melting_rate,'type','node','layer',1); 
			self.perturbation_melting_rate=project3d(md,'vector',self.perturbation_melting_rate,'type','node','layer',1); 
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','node','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end

			if isnan(self.floatingice_melting_rate),
				self.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.floatingice_melting_rate specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{
			%Set default parameter values
			self.Kh = 25; % unit: m2 s-1
			self.Ah = 25; % unit: m2 s-1
			self.Dmin   = 2.8; % unit: m
			self.Utide  = 0.01; % unit: m s-1
			self.Cd= 2.5e-3; % unit: -
			self.Cd_top=1.1e-3; % unit -
			self.f_cori = 1.37e-4; % unit: s-1

			%Use entrainment with Holland et al. (2006).
			self.isentrainment = 0;
			self.ismelt = 1;

			%Time stepping specials
			self.subtimestep = 72; % unit: s-1
			self.diagnostic_frequency=1; % unit: -

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.ismasstransport==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			end
			if ismember('ThermalAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.isthermal==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			end

			%Check consistency for initial state of D, T, S, vx, vy
			md = checkfield(md,'fieldname','basalforcings.D','size','NaN',1,'Inf',1,'size','[d.mesh.numberofvertices, 1]);
			md = checkfield(md,'fieldname','basalforcings.vx','size','NaN',1,'Inf',1,'size','[d.mesh.numberofvertices, 1]);
			md = checkfield(md,'fieldname','basalforcings.vy','size','NaN',1,'Inf',1,'size','[d.mesh.numberofvertices, 1]);
			md = checkfield(md,'fieldname','basalforcings.T','size','NaN',1,'Inf',1,'size','[d.mesh.numberofvertices, 1]);
			md = checkfield(md,'fieldname','basalforcings.S','size','NaN',1,'Inf',1,'size','[d.mesh.numberofvertices, 1]);

			%Check consistency for forcing values
			md = checkfield(md,'fieldname','basalforcings.forcing_depth','NaN',1,'Inf',1,'size',[1,NaN],'<=',0);
			md = checkfield(md,'fieldname','basalforcings.forcing_temperature','size',[1,1,numel(md.basalforcings.forcing_depth)]);
			md = checkfield(md,'fieldname','basalforcings.forcing_salinity','size',[1,1,numel(md.basalforcings.forcing_depth)]);
			for i=1:numel(md.basalforcings.forcing_depth)
				md = checkfield(md,'fieldname',['basalforcings.forcing_temperature{' num2str(i) '}'],'field',md.basalforcings.forcing_temperature{i},'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'>=',0,'timeseries',1);
				md = checkfield(md,'fieldname',['basalforcings.forcing_salinity{' num2str(i) '}'],'field',md.basalforcings.forcing_salinity{i},'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'>=',0,'timeseries',1);
			end

			%Check consistency for parameters
			md = checkfield(md,'fieldname','Kh','numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','Ah','numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','Dmin','numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','Utide','numel',1,'NaN',1,'Inf',1);

			md = checkfield(md,'fieldname','isentrainment','numel',1,'NaN',1,'Inf',1,'values',[0,1]);
			md = checkfield(md,'fieldname','ismelt','numel',1,'NaN',1,'Inf',1,'values',[0,1]);
			md = checkfield(md,'fieldname','Cd','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','Cd_top','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','f_cori','numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','subtimestep','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','diagnostic_frequency','numel',1,'NaN',1,'Inf',1,'>',0);
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',10,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts)
			WriteData(fid,prefix,'object',self,'fieldname','D','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','vx','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','vy','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','temperature','format','DoubleMat','mattype',1)
			WriteData(fid,prefix,'object',self,'fieldname','salinity','format','DoubleMat','mattype',1)

			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);

			%Write forcing values
			WriteData(fid,prefix,'object',self,'fieldname','forcing_depth','format','DoubleMat','name','md.basalforcings.forcing_depth');
			WriteData(fid,prefix,'object',self,'fieldname','forcing_temperature','format','MatArray','name','md.basalforcings.forcing_temperature','timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','forcing_salinity','format','MatArray','name','md.basalforcings.forcing_salinity','timeserieslength',md.mesh.numberofvertices+1,'yts',yts);

			%Write parameters
			WriteData(fid,prefix,'object',self,'fieldname','Ah','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','Kh','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','Dmin','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','Utide','format','Double');

			WriteData(fid,prefix,'object',self,'fieldname','isentrainment','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','ismelt','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','ismelt','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','Cd','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','Cd_top','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','f_cori','format','Double');

			WriteData(fid,prefix,'object',self,'fieldname','subtimestep','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','diagnostic_frequency','format','Integer');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.basalforcings.groundedice_melting_rate'],self.groundedice_melting_rate);
			writejs1Darray(fid,[modelname '.basalforcings.floatingice_melting_rate'],self.floatingice_melting_rate);
			writejs1Darray(fid,[modelname '.basalforcings.perturbation_melting_rate'],self.perturbation_melting_rate);
			writejs1Darray(fid,[modelname '.basalforcings.geothermalflux'],self.geothermalflux);

		end % }}}
	end
end
