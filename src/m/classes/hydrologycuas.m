%HYDROLOGYCUAS class definition
%
%   Usage:
%      hydrologycuas=hydrologycuas();

classdef hydrologycuas
	properties (SetAccess=public)
		head = NaN;
		transmissivity = NaN;
		conductivity = NaN;
		layer_thickness = NaN;

		bump_height = NaN;
		bump_spacing= NaN;

		Tmin  = NaN;
		Tmax  = NaN;
		
		ss = NaN;
		sy = NaN;
		englacial_input = NaN;
		moulin_input = NaN;

		% Parameters for channels
		ischannel_creep =0;
		ischannel_melt  =0;
		ischannel_cavity=0;
		isconfined=0;
		melt_flag=0;
		unconfinedSmooth=0;

		% Boundary conditions
		spchead = NaN;
		neumannflux = NaN;

		steps_per_step=NaN;
		averaging=NaN;

		requested_outputs = {};
	end
	methods
		function self = extrude(self,md) % {{{
			self.head=project3d(md,'vector',self.spcsediment_head,'type','node','layer',1);
		end % }}}
		function self = hydrologycuas(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(self,varargin{1});
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end% }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'HydrologyHead','HydrologyTransmissivityEffective','HydrologyStorage','HydrologyTransmissivity',...
				'EffectivePressure'};
		end % }}}
		function self = initialize(self,md) % {{{
			nV = md.mesh.numberofvertices;

			self.layer_thickness= 1*ones(nV,1);
			self.conductivity   = 0.001*ones(nV,1);
			self.transmissivity = self.transmissivity.*self.layer_thickness;

			self.ss             = 0.0000982977696*ones(nV,1);
			self.sy             = 0.4*ones(nV,1);

			% Source term...
			self.englacial_input = zeros(nV,1);
			self.moulin_input    = zeros(nV,1);

			% Opening factor.
			self.bump_height  = 0.1*ones(nV,1);
			self.bump_spacing = 10.0*ones(nV,1);

			self=setdefaultparameters(self);
			return;
		end % }}}
		function self = setdefaultparameters(self)% {{{
			%Parameters from original CUAS document
			self.Tmax           = 20.0;
			self.Tmin           = 1e-7;
			
			%self.ss        = 0.0000982977696;
			%self.sy        = 0.4;

			self.ischannel_creep = 1;
			self.ischannel_melt  = 1;
			self.ischannel_cavity= 1;
			self.isconfined= 1;

			self.melt_flag=1;

			self.unconfinedSmooth=0.0;

			self.requested_outputs        = {'default'};
		end  % }}}
		function md = checkconsistency(self,md,solution,analyses)% {{{
			%Early return
			if ~ismember('HydrologyCuasAnalysis',analyses) & ~ismember('HydrologyCuasAnalysis',analyses),
				return;
			end
			md = checkfield(md,'fieldname','hydrology.head','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.conductivity','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.layer_thickness','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.transmissivity','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);

			md = checkfield(md,'fieldname','hydrology.Tmin','numel',1,'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.Tmax','numel',1,'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.ss','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.sy','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);

			md = checkfield(md,'fieldname','hydrology.moulin_input','NaN',1,'Inf',1,'timeseries',1);

			md = checkfield(md,'fieldname','hydrology.neumannflux','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.spchead','Inf',1,'timeseries',1);

			md = checkfield(md,'fieldname','hydrology.ischannel_creep','values',[0,1],'Inf',1,'NaN',1,'numel',1);
			md = checkfield(md,'fieldname','hydrology.ischannel_melt','values',[0,1],'Inf',1,'NaN',1,'numel',1);
			md = checkfield(md,'fieldname','hydrology.ischannel_cavity','values',[0,1],'Inf',1,'NaN',1,'numel',1);
			md = checkfield(md,'fieldname','hydrology.isconfined','values',[0,1],'Inf',1,'NaN',1,'numel',1);
		end% }}}
		function disp(self)% {{{
			disp(sprintf('   hydrology CUAS (Confined-Unconfined Aquifer System) parameters:'));
			fielddisplay(self,'head','initial hydraulic head (m)');
			fielddisplay(self,'conductivity','Conductivity of layer (m s-1)');
			fielddisplay(self,'layer_thickness','EPM layer thickness (m)');
			fielddisplay(self,'transmissivity','transmissivity of EPM layer (m^2 s-1)');
			fielddisplay(self,'Tmin','Minimum transmissivity to be allowed in the evolution (m^2 s-1)');
			fielddisplay(self,'Tmax','Maximum transmissivity to be allowed in the evolution (m^2 s-1)');
			fielddisplay(self,'ss','specific storage, ss (unit: m^-1)');
			fielddisplay(self,'sy','specific yield, Sy (unit: 1)');
			fielddisplay(self,'moulin_input','liquid water input from moulins (at the vertices) to subglacial system (m^3/s)');

			fielddisplay(self,'ischannel_creep','Evolve channel due to creep (default: 1)');
			fielddisplay(self,'ischannel_melt','Evolve channel due to melting (default: 1)');
			fielddisplay(self,'ischannel_cavity','Evolve channel due to cavity (default: 1)');
			fielddisplay(self,'isconfined','Select CUAS model to solve confined or unconfined. 0: unconfined, 1: confined. (default: 1)');
			fielddisplay(self,'unconfinedSmooth','Smoothing factor for unconfined aquifer system (default: 0.0)');
		end % }}}
		function marshall(self,prefix,md,fid)% {{{
			WriteData(fid,prefix,'name','md.hydrology.model','data',8,'format','Integer');
			yts = md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','head','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','conductivity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','layer_thickness','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','transmissivity','format','DoubleMat','mattype',1);
			
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','Tmin','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','Tmax','format','Double');

			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ss','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sy','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','englacial_input','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','moulin_input','format','DoubleMat','mattype',1);

			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_height','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_spacing','format','DoubleMat','mattype',1);


			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ischannel_cavity','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ischannel_melt','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ischannel_creep','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','isconfined','format','Boolean');

			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','melt_flag','format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','unconfinedSmooth','format','Double');

			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','spchead','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','neumannflux','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);

			% requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];  %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end% }}}
	end
end
