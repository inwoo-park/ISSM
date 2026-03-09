%HYDROLOGYCUAS class definition
%
%   Usage:
%      hydrologycuas=hydrologycuas();

classdef hydrologycuas
	properties (SetAccess=public)
		conductivity = NaN;
		Tmin  = NaN;
		Tmax  = NaN;
		Tinit = NaN;

		ischannel_creep =0;
		ischannel_melt  =0;
		ischannel_cavity=0;
		isconfined=0;

		steps_per_step=NaN;
		averaging=NaN;
	end
	methods
		function self = extrude(self,md)    % {{{
			self.head=project3d(md,'vector',self.spcsediment_head,'type','node','layer',1);
		end    % }}}
		function self = hydrologycuas(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end% }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SedimentHead','SedimentHeadResidual','EffectivePressure'};
			if self.isefficientlayer,
				list=[list,{'EplHead','HydrologydcMaskEplactiveNode','HydrologydcMaskEplactiveElt','EplHeadSlopeX','EplHeadSlopeY','HydrologydcEplThickness'}];
			end
			if self.steps_per_step>1 | self.step_adapt,
				list = [list,'EffectivePressureSubstep','SedimentHeadSubstep'];
				if self.isefficientlayer,
					list = [list,'EplHeadSubstep','HydrologydcEplThicknessSubstep'];
				end
			end
		end % }}}
		function self = initialize(self,md) % {{{
			self.epl_colapse_thickness = self.sediment_transmitivity/self.epl_conductivity;
			if isnan(self.basal_moulin_input),
				self.basal_moulin_input=zeros(md.mesh.numberofvertices,1);
				disp('      no hydrology.basal_moulin_input specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self)% {{{
			%Parameters from ?? (original CUAS document)
			self.Tmin      = 1e-7;
			self.Tmax      = 20.0;
			self.Tinit     = 0.2;
			
			self.ss        = 0.0000982977696;
			self.sy        = 0.4;

			self.ischannel_creep = 1;
			self.ischannel_melt  = 1;
			self.ischannel_cavity= 1;
			self.isconfined= 1;

			% Opening factor.
			self.roughness = 1.0;

			self.requested_outputs        = {'default'};
		end  % }}}
		function md = checkconsistency(self,md,solution,analyses)% {{{
			%Early return
			if ~ismember('HydrologyCuasAnalysis',analyses) & ~ismember('HydrologyCuasAnalysis',analyses),
				return;
			end
			md = checkfield(md,'fieldname','hydrology.conductivity','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.Tmin','numel',1,'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.Tmax','numel',1,'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.Tinit','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.ss','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.sy','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.roughness','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1,'>=',0);

			md = checkfield(md,'fieldname','hydrology.ischannel_creep','values',[0,1],'Inf',1,'NaN',1,'numel',1);
			md = checkfield(md,'fieldname','hydrology.ischannel_melt','values',[0,1],'Inf',1,'NaN',1,'numel',1);
			md = checkfield(md,'fieldname','hydrology.ischannel_cavitiy','values',[0,1],'Inf',1,'NaN',1,'numel',1);
			md = checkfield(md,'fieldname','hydrology.isconfined','values',[0,1],'Inf',1,'NaN',1,'numel',1);
		end% }}}
		function disp(self)% {{{
			disp(sprintf('   hydrology CUAS (Confined-Unconfined Aquifer System) parameters:'));
			fielddisplay(self,'conductivity','Conductivity of layer (m s-1)');
			fielddisplay(self,'Tmin','Minimum transmissivity to be allowed in the evolution (m^2 s-1)');
			fielddisplay(self,'Tmax','Maximum transmissivity to be allowed in the evolution (m^2 s-1)');
			fielddisplay(self,'Tinit','Initial transmissivity to be used in the evolution (m^2 s-1)');
			fielddisplay(self,'ss','specific storage, ss (unit: m^-1)');
			fielddisplay(self,'sy','specific yield, Sy (unit: 1)');
			fielddisplay(self,'roughness','roughness factor for opening term (default: 1.0)');

			fielddisplay(self,'ischannel_creep','Evolve channel due to creep (default: 1)');
			fielddisplay(self,'ischannel_melt','Evolve channel due to melting (default: 1)');
			fielddisplay(self,'ischannel_cavitiy','Evolve channel due to cavity (default: 1)');
			fielddisplay(self,'isconfined','Select CUAS model to solve confined or unconfined. 0: unconfined, 1: confined. (default: 1)');
		end % }}}
		function marshall(self,prefix,md,fid)% {{{
			WriteData(fid,prefix,'name','md.hydrology.model','data',8,'format','Integer');

			WriteData(fid,prefix,'object',self,'fieldname','head','format','Double');



			WriteData(fid,prefix,'object',self,'fieldname','water_compressibility','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','isefficientlayer','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','penalty_factor','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','penalty_lock','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','rel_tol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','max_iter','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','step_adapt','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','sedimentlimit_flag','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','transfer_flag','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','unconfined_flag','format','Integer');
			if self.sedimentlimit_flag==1,
				WriteData(fid,prefix,'object',self,'fieldname','sedimentlimit','format','Double');
			end
			if self.transfer_flag==1,
				WriteData(fid,prefix,'object',self,'fieldname','leakage_factor','format','Double');
			end
			WriteData(fid,prefix,'object',self,'fieldname','basal_moulin_input','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)

			WriteData(fid,prefix,'object',self,'fieldname','spcsediment_head','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','sediment_compressibility','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','sediment_porosity','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','sediment_thickness','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','sediment_transmitivity','format','DoubleMat','mattype',1');
			WriteData(fid,prefix,'object',self,'fieldname','mask_thawed_node','format','DoubleMat','mattype',1);
			if self.isefficientlayer==1,
				WriteData(fid,prefix,'object',self,'fieldname','spcepl_head','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'fieldname','mask_eplactive_node','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'fieldname','epl_compressibility','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_porosity','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_initial_thickness','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_colapse_thickness','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_thick_comp','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','epl_max_thickness','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_conductivity','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','eplflip_lock','format','Integer');
			end
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
