%DAMAGEICE class definition
%
%   Usage:
%      damage=damage();

classdef damage
	properties (SetAccess=public)  
		%damage 
		isdamage            = 0;
		D                   = NaN;
		law                 = 0;
		spcdamage           = NaN; 
		max_damage          = 0;
	
		%numerical
		stabilization       = 0;
		maxiter             = 0;
		elementinterp       = '';
		
		%general parameters for evolution law: 
		stress_threshold    = 0;
		stress_ubound       = 0;
		kappa               = 0;
		c1                  = 0;
		c2                  = 0;
		c3                  = 0;
		c4                  = 0;
		healing             = 0;
		isdamage_exponent   = 0;
		ispressure_ssa      = 0;
		isPeff              = 0;
		equiv_stress		  = 0;
		requested_outputs   = {};

		equiv_stress_alpha  = 0; % For hayhurst
		equiv_stress_beta   = 0; % For hayhurst
		equiv_stress_mu     = 0; % For Coulomb
	end
	methods
		function self = damage(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('damage');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2)
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.D=project3d(md,'vector',self.D,'type','node');
			self.spcdamage=project3d(md,'vector',self.spcdamage,'type','node');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%damage parameters: 
			self.isdamage=0;
			self.D=0;
			self.law=0;
			
			self.max_damage=1-1e-5; %if damage reaches 1, solve becomes singular, as viscosity becomes nil
		
			%Type of stabilization used
			self.stabilization=4;
			
			%Maximum number of iterations
			self.maxiter=100;

			%finite element interpolation
			self.elementinterp='P1';

			%damage evolution parameters 
			self.stress_threshold=1.3e5;
			self.kappa=2.8;
			self.healing=0;
			self.c1=0;
			self.c2=0;
			self.c3=0;
			self.c4=0;
			self.equiv_stress=0;

			self.isdamage_exponent=0;
			self.ispressure_ssa=0;
			self.isPeff=0;

			% Criterion of Hayhurst (Pralong et al., 2005)
			self.equiv_stress_alpha=0.21;
			self.equiv_stress_beta=0.63;

			% Criterion of Coulomb (Vaughan 1993)
			self.equiv_stress_mu = 0.1;

			%output default:
			self.requested_outputs={'default'};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			
			md = checkfield(md,'fieldname','damage.isdamage','values',[1,0]);
			if self.isdamage
				%Check class of md.materials.
				if ~isa(md.materials,'matdamageice')
					error('Error: Invalid class of md.materials (=%s). md.materials should be a subclass of "matdamageice".',class(md.materials));
				end

				md = checkfield(md,'fieldname','damage.law','numel',[1],'values',[0,1,2,3,4,5]);
				md = checkfield(md,'fieldname','damage.D','>=',0,'<=',self.max_damage,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','damage.spcdamage','Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','damage.max_damage','<',1,'>=',0);
				md = checkfield(md,'fieldname','damage.stabilization','numel',[1],'values',[0 1 2 4]);
				md = checkfield(md,'fieldname','damage.maxiter','>=',0);
				md = checkfield(md,'fieldname','damage.elementinterp','values',{'P1','P2'});
				md = checkfield(md,'fieldname','damage.stress_threshold','>=',0);
				md = checkfield(md,'fieldname','damage.stress_ubound','>=',0);
				md = checkfield(md,'fieldname','damage.kappa','>',1);
				md = checkfield(md,'fieldname','damage.healing','>=',0);
				md = checkfield(md,'fieldname','damage.c1','>=',0);
				md = checkfield(md,'fieldname','damage.c2','>=',0);
				md = checkfield(md,'fieldname','damage.c3','>=',0);
				md = checkfield(md,'fieldname','damage.c4','>=',0);
				md = checkfield(md,'fieldname','damage.isdamage_exponent','numel',[1],'values',[0 1 2]);
				md = checkfield(md,'fieldname','damage.ispressure_ssa','numel',[1],'values',[0 1 2]);
				md = checkfield(md,'fieldname','damage.isPeff','numel',[1],'values',[0 1]);
				md = checkfield(md,'fieldname','damage.equiv_stress','numel',[1],'values',[0 1 2 3]);
				md = checkfield(md,'fieldname','damage.requested_outputs','stringrow',1);
			elseif (self.law~=0)
				if (strcmp(solution,'DamageEvolutionSolution'))
					error('Invalid evolution law (md.damage.law) for a damage solution');
				end
			end
		end % }}}
		function list=defaultoutputs(self,md) % {{{

			if strcmp(domaintype(md.mesh),'2Dhorizontal')
            list = {'DamageDbar'};
         else
            list = {'DamageD'};
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Damage:\n'));

			fielddisplay(self,'isdamage','is damage mechanics being used? {true,false}');
			if self.isdamage
				fielddisplay(self,'law','damage law {''0: analytical'',''1: pralong'',''4: test'',''5: linear''}');
				fielddisplay(self,'D','damage tensor (scalar)');
				fielddisplay(self,'spcdamage','damage constraints (NaN means no constraint)');
				fielddisplay(self,'max_damage','maximum possible damage (0<=max_damage<1)');
				
				fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: SUPG (not working), 4: flux corrected transport');
				fielddisplay(self,'maxiter','maximum number of non linear iterations');
				fielddisplay(self,'elementinterp','interpolation scheme for finite elements {''P1'',''P2''}');
				fielddisplay(self,'stress_threshold','stress threshold for damage initiation (Pa)');
				fielddisplay(self,'stress_ubound','stress upper bound for damage healing (Pa), arctan law');
				fielddisplay(self,'kappa','ductility parameter for stress softening and damage');
				if self.law < 4
					fielddisplay(self,'c1','damage parameter 1');
					fielddisplay(self,'c2','damage parameter 2');
					fielddisplay(self,'c3','damage parameter 3');
					fielddisplay(self,'c4','damage parameter 4');
				elseif self.law == 4
					fprintf('\n');
					disp(sprintf('   Damage source term in Pralong et al. (2005)'));
					disp(sprintf('   f_D = c1 * max((sigma_cr - sigma_th),0)^c2 (1-D)^{k_sigma}'));
					fprintf('\n');
					disp(sprintf('   exponent in damage term, k_sigma'));
					disp(sprintf('   0 : c3'));
					disp(sprintf('   1 : k1 * (sqrt(max(s_inv1,0)) - healing * sqrt(max(-s_inv1,0)))'));
					disp(sprintf('   2 : k1 + k2 * s_inv1'));
					fielddisplay(self,'c1','damage parameter 1');
					fielddisplay(self,'c2','damage parameter 2');
					fielddisplay(self,'c3','damage parameter 3');
					fielddisplay(self,'isdamage_exponent','damage exponent parameter (k_sigma). 0: constant, 1: Pralong (2005), 2: Duddu et al. (2020)');
					fielddisplay(self,'ispressure_ssa','pressure for SSA2D. 0: surface, 1: mid, 2: bottom');
				elseif slef.law == 5
					fielddisplay(self,'c1','damage parameter 1 ( F = c1 * max( stress_equiv - stress_threshold, 0)');
				end
				fielddisplay(self,'healing','damage healing parameter');
				fielddisplay(self,'equiv_stress','0: von Mises, 1: max prinecipal, 2: Hayhurst criterion 3: Coulomb');
				if self.equiv_stress == 2
					disp(sprintf('\n   Hayhurst criterion'));
					disp(sprintf('      sigma_equiv = alpha*sigma_1 + beta*sigma_{VM} + (1-alpha-beta)*(sigma_1+sigma_2+sigma_3)'));
					disp(sprintf('      sigma_i : principal stress value'));
					fielddisplay(self,'equiv_stress_alpha','alpha parameter for Hayhurst criterion (default: 0.21)');
					fielddisplay(self,'equiv_stress_beta','beta parameter for Hayhurst criterion (default: 0.63)');
				elseif self.equiv_stress == 3
					disp(sprintf('\n   Coulomb criterion'));
					disp(sprintf('      sigma_equiv = (sigma_1 - sigma_3)/2 + mu (sigma_1 + sigma_3)/2'));
					disp(sprintf('      sigma_i : principal stress value'));
					fielddisplay(self,'equiv_stress_mu','mu parameter for Coulomb criterion (default: 0.1)');
				end
				fprintf('\n');
				fielddisplay(self,'isPeff','considering water pressure at base. 0: off, 1: on (default: 0)');
				fielddisplay(self,'requested_outputs','additional outputs requested');
			end

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
		
			WriteData(fid,prefix,'object',self,'fieldname','isdamage','format','Boolean');
			if self.isdamage
				WriteData(fid,prefix,'object',self,'fieldname','law','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','D','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'fieldname','spcdamage','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'fieldname','max_damage','format','Double');

				WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','maxiter','format','Integer');
				WriteData(fid,prefix,'name','md.damage.elementinterp','data',self.elementinterp,'format','String');
				WriteData(fid,prefix,'object',self,'fieldname','stress_threshold','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','stress_ubound','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','kappa','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','c1','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','c2','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','c3','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','c4','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','healing','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','equiv_stress','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','equiv_stress_alpha','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','equiv_stress_beta','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','equiv_stress_mu','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','isdamage_exponent','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','ispressure_ssa','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','isPeff','format','Boolean');
			end

			%process requested outputs
			outputs = self.requested_outputs;
			pos = find(ismember(outputs,'default'));
			if ~isempty(pos)
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			if self.isdamage
				WriteData(fid,prefix,'data',outputs,'name','md.damage.requested_outputs','format','StringArray');
			end

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.damage.isdamage'],self.isdamage);
			if self.isdamage
				error('savemodeljs error message: not implemented  yet!');
			end

		end % }}}
	end
end
