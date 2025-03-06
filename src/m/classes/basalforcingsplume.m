%BASAL FORCINGS class definition
%
%   Usage:
%      basalforcings=basalforcingsplume();

classdef basalforcingsplume
	properties (SetAccess=public) 
		groundedice_melting_rate  = NaN;
		geothermalflux            = NaN;
		D = NaN;
		vx= NaN;
		vy= NaN;
		temperature = NaN;
		salinity = NaN;

		% Dirichlet boundary
		spcD  = NaN;
		spcvx = NaN;
		spcvy = NaN;
		spctemperature = NaN;
		spcsalinity    = NaN;

		% Neumann boundary

		% parameters
		entrainment = 0;
	end
	methods
		function self = basalforcingsplume(varargin) % {{{
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
			fielddisplay(self,'vx','initial plume velocity in x-direction [m s-1].');
			fielddisplay(self,'vy','initial plume velocity in y-direction [m s-1].');
			fielddisplay(self,'temperature','initial plume temperature [K].');
			fielddisplay(self,'salinity','initial plume salinity [psu].');

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
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',10,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts)
			WriteData(fid,prefix,'object',self,'fieldname','D','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts)
			WriteData(fid,prefix,'object',self,'fieldname','vx','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts)
			WriteData(fid,prefix,'object',self,'fieldname','vy','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts)
			WriteData(fid,prefix,'object',self,'fieldname','temperature','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts)
			WriteData(fid,prefix,'object',self,'fieldname','salinity','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts)

			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.basalforcings.groundedice_melting_rate'],self.groundedice_melting_rate);
			writejs1Darray(fid,[modelname '.basalforcings.floatingice_melting_rate'],self.floatingice_melting_rate);
			writejs1Darray(fid,[modelname '.basalforcings.perturbation_melting_rate'],self.perturbation_melting_rate);
			writejs1Darray(fid,[modelname '.basalforcings.geothermalflux'],self.geothermalflux);

		end % }}}
	end
end
