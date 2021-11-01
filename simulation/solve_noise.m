%% set different kinds of miRNA regulations
% This code can simulate miRNA regulations with
% (1) competing RNAs
% (2) repetitive targets of same miRNAs
% (3) multiple targets of different miRNAs
% simulation type
type = 1;

%% set fixed parameters
% parameters should be set according to the study purpose
% for all situations (No. 1-10)
kR = 5e-3;      gR = 1e-4;      kT = 1e-3;      gT = gR;
kon = 1e-4;     koff = 5e-5;    alpha = 0.5;    gC = 2 * gR;
kp = 1e-2;      gp = 5e-6;

% for simulation type 1 (No. 11-16)
kT2 = kT;      gT2 = gR;       kon2 = kon;     koff2 = koff; 
alpha2 = 0.5;  gC2 = 2 * gR;

% for simulation type 2 (No. 17-18)
konDouble = kon;     koffDouble = koff;

% for simulation type 3 (No. 19-27)
kR2 = kR;       gR2 = gR;
konB = kon;     koffB = koff;
konA_CB = kon;  koffA_CB = koff;
konB_CA = kon;  koffB_CA = koff;
gCB = gC;       alphaB = alpha;

% for simulation type 2 & 3 (No. 28)
gCAB = gC;      alphaAB = alpha;

%% set altered parameters
% here we take the study of strong and weak competing RNAs as an example
% assign which parameters will be altered according to the parameter No.
% The first one should be the control condition

% samples for type == 1
alter_para_No = [1 11 14]; % kR kT2 koff2

alter_para = [  kR * [0 1 2 2]; 
                kT2 * [0 0 10 15];
                koff2 * [1 1 1 100]];

% samples for type == 2
%{
alter_para_No = [1 5 17]; % kR kon konDouble

alter_para = [  kR * [0 1 1 1];
                kon * [1 1 1 2];
                konDouble * [0 0 1 0]];
%}

% samples for type == 3
%{
alter_para_No = [1 19]; % kR kR2

alter_para = [  kR * [0 1 1 2];
                kR2 * [0 0 1 0]];
%}
            
%% set simulation parameters
% the number of different expression levels of the observed gene
N_level = 100;
% the range of log10(kT1)
min_KT1 = -3;
max_kT1 = -1;

%% start to run
% create kT list
kT_list = 10.^linspace(min_KT1, max_kT1, N_level);

% the number of different simulation conditions
[N_para, N_alter] = size(alter_para);

% create parameter matrix
all_para = [kR gR kT gT kon koff alpha gC kp gp ...
    kT2 gT2 kon2 koff2 alpha2 gC2 ...
    konDouble koffDouble kR2 gR2 ...
    konB koffB konA_CB koffA_CB konB_CA koffB_CA ...
    gCB alphaB gCAB alphaAB];
para_mat = repmat(all_para, 4, 1)';
for i = 1:N_para
    para_mat(alter_para_No(i), :) = alter_para(i, :);
end
    
% create result matrix
exp_mean = zeros(N_alter, N_level);
exp_CV = zeros(N_alter, N_level);

% start
for alter_i = 1:1:N_alter
    this_para = para_mat(:, alter_i);
    
    for level_i = 1:1:N_level
        this_para(3) = kT_list(level_i);
        % solve steady state
        INIT = steady_state(type, this_para);
        INIT = num2cell(INIT);
        % the steady state vector N: R T C(CA) p T2 C2 CB CAB R2
        [R, T, C, p, T2, C2, CB, CAB, R2] = deal(INIT{:});
        
        this_para2 = num2cell(this_para);
        [KR, GR, KT, GT, Kon, Koff, Alpha, GC, Kp, Gp, ...
         KT2, GT2, Kon2, Koff2, Alpha2, GC2, ...
         KonDouble, KoffDouble, KR2, GR2, ...
         KonB, KoffB, KonA_CB, KoffA_CB, KonB_CA, KoffB_CA, ...
         GCB, AlphaB, GCAB, AlphaAB] = deal(this_para2{:});
                
        exp_mean(alter_i, level_i) = p;
          
        % transition rate vector
        vTp = [KT, GT*T, Kp*T, Gp*p];
        vT2 = [KT2, GT2*T2];
        vR = [KR, GR*R];
        vR2 = [KR2, GR2*R2];
        vC = [Kon*T*R, Koff*C, GC*(1 - Alpha)*C, GC*Alpha*C];
        vC2 = [Kon2*T2*R, Koff2*C2, GC2*(1 - Alpha2)*C2, GC2*Alpha2*C2];
        vCB_rep = [Kon*T*R, Koff*CB, GCB*(1 - AlphaB)*CB, GCB*AlphaB*CB];
        vCB_com = [KonB*T*R2, KoffB*CB, GCB*(1 - AlphaB)*CB, GCB*AlphaB*CB];
        vCAB_rep = [KonDouble*C*R, KonDouble*CB*R, KoffDouble*CAB, KoffDouble*CAB, GCAB*(1 - AlphaAB)*CAB, GCAB*AlphaAB*CAB];
        vCAB_com = [KonB_CA*C*R2, KonA_CB*CB*R, KoffB_CA*CAB, KoffA_CB*CAB, GCAB*(1 - AlphaAB)*CAB, GCAB*AlphaAB*CAB];
        
        % stoichiometrix matrix 
            %[  R   T   C   p   T2  C2  CB  CAB R2]
        STp = [	0	1	0	0	0	0	0	0	0;
            	0   -1	0	0	0	0	0	0	0;
                0   0	0	1	0	0	0	0	0;
                0   0	0	-1	0	0	0	0	0];

        ST2 = [	0	0	0	0	1	0	0	0	0;
            	0   0	0	0	-1	0	0	0	0];
            
        SR = [	1	0	0	0	0	0	0	0	0;
            	-1   0	0	0	0	0	0	0	0];

        SR2 = [	0	0	0	0	0	0	0	0	1;
            	0   0	0	0	0	0	0	0	-1];
            
        SC = [	-1	-1	1	0	0	0	0	0	0;
            	1   1	-1	0	0	0	0	0	0;
                1	0	-1	0	0	0	0	0	0;
            	0   0	-1	0	0	0	0	0	0];     
 
        SC2 = [	-1	0	0	0	-1	1	0	0	0;
            	1   0	0	0	1	-1	0	0	0;
                1	0	0	0	0	-1	0	0	0;
            	0   0	0	0	0	-1	0	0	0];
            
        SCB_rep = [	
                -1	-1	0	0	0	0	1	0	0;
            	1   1	0	0	0	0	-1	0	0;
                1	0	0	0	0	0	-1	0	0;
            	0   0	0	0	0	0	-1	0	0]; 
            
        SCB_com = [	
                0	-1	0	0	0	0	1	0	-1;
            	0   1	0	0	0	0	-1	0	1;
                1	0	0	0	0	0	-1	0	0;
            	0   0	0	0	0	0	-1	0	0];   
            
        SCAB_rep = [
                -1	0	-1	0	0	0	0	1	0;
            	-1  0	0	0	0	0	-1	1	0;
                1	0	1	0	0	0	0	-1	0;
            	1   0	0	0	0	0	1	-1	0;
                2   0	0	0	0	0	0	-1	0;
                0   0	0	0	0	0	0	-1	0];  
            
        SCAB_com = [
                0	0	-1	0	0	0	0	1	-1;
            	-1  0	0	0	0	0	-1	1	0;
                0	0	1	0	0	0	0	-1	1;
            	1   0	0	0	0	0	1	-1	0;
                1   0	0	0	0	0	0	-1	1;
                0   0	0	0	0	0	0	-1	0];  
            
   
        % dv = dv/dN
            %  [    R   T   C   p   T2  C2  CB  CAB R2]   
        dvTp = [	0	0	0	0	0	0	0	0	0;
                    0   GT	0	0	0	0	0	0	0;
                    0   Kp	0	0	0	0	0	0	0;
                    0   0	0	Gp	0	0	0	0	0];
 
        dvT2 = [    0	0	0	0	0	0	0	0	0;
                    0   0	0	0	GT2	0	0	0	0];
            
        dvR = [     0	0	0	0	0	0	0	0	0;
                    GR  0	0	0	0	0	0	0	0];

        dvR2 = [	0	0	0	0	0	0	0	0	0;
                    0   0	0	0	0	0	0	0	GR2];
                
        dvC = [     Kon*T	Kon*R	0	0	0	0	0	0	0;
                    0   0	Koff	0	0	0	0	0	0;
                    0	0	GC*(1 - Alpha)	0	0	0	0	0	0;
                    0   0	GC*Alpha	0	0	0	0	0	0];     
 
        dvC2 = [	Kon2*T2	0	0	0	Kon2*R	0	0	0	0;
                    0   0	0	0	0	Koff2	0	0	0;
                    0	0	0	0	0	GC2*(1 - Alpha2)	0	0	0;
                    0   0	0	0	0	GC2*Alpha2	0	0	0];
            
        dvCB_rep = [	
                    Kon*T	Kon*R	0	0	0	0	0	0	0;
                    0   0	0	0	0	0	Koff	0	0;
                    0	0	0	0	0	0	GCB*(1 - AlphaB)	0	0;
                    0   0	0	0	0	0	GCB*AlphaB	0	0]; 
  
        dvCB_com = [	
                    0	KonB*R2	0	0	0	0	0	0	KonB*T;
                    0   0	0	0	0	0	KoffB	0	0;
                    0	0	0	0	0	0	GCB*(1 - AlphaB)	0	0;
                    0   0	0	0	0	0	GCB*AlphaB	0	0];
                
                
            %[  R   T   C   p   T2  C2  CB  CAB R2]
        dvCAB_rep = [
                    KonDouble*C	0	KonDouble*R	0	0	0	0	0	0;
                    KonDouble*CB  0	0	0	0	0	KonDouble*R	0	0;
                    0	0	0	0	0	0	0	KoffDouble	0;
                    0   0	0	0	0	0	0	KoffDouble	0;
                    0   0	0	0	0	0	0	GCAB*(1 - AlphaAB)	0;
                    0   0	0	0	0	0	0	GCAB*AlphaAB	0];  
            
        dvCAB_com = [
                    0	0	KonB_CA*R2	0	0	0	0	0	KonB_CA*C;
                    KonA_CB*CB  0	0	0	0	0	KonA_CB*R	0	0;
                    0	0	0	0	0	0	0	KoffB_CA	0;
                    0   0	0	0	0	0	0	KoffA_CB	0;
                    0   0	0	0	0	0	0	GCAB*(1 - AlphaAB)	0;
                    0   0	0	0	0	0	0	GCAB*AlphaAB	0];
        
        if type == 1
            v = [vTp vT2 vR vC vC2]';
            S = [STp; ST2; SR; SC; SC2]';
            dv = [dvTp; dvT2; dvR; dvC; dvC2]';
            S = S(1:6, :);
            dv = dv(1:6, :);
        elseif type == 2
            v = [vTp vR vC vCB_rep vCAB_rep]';
            S = [STp; SR; SC; SCB_rep; SCAB_rep]';
            dv = [dvTp; dvR; dvC; dvCB_rep; dvCAB_rep]';
            S = S([1:4,7:8], :);
            dv = dv([1:4,7:8], :);
        else
            v = [vTp vR vR2 vC vCB_com vCAB_com]';
            S = [STp; SR; SR2; SC; SCB_rep; SCAB_com]';
            dv = [dvTp; dvR; dvR2; dvC; dvCB_rep; dvCAB_com]';
            S = S([1:4,7:9], :);
            dv = dv([1:4,7:9], :);
        end
        
    
        N_reaction = length(v);
        if length(size(exp_CV)) == 2
            exp_CV = zeros(N_alter, N_level, N_reaction);
        end
        % divide the noise by different reactions
        % when reaction = 0, calculate the total noise
        for reaction_i = 0:1:N_reaction
            if reaction_i == 0
                v1 = v;
            else
                v1 = zeros(size(v));
                v1(reaction_i) = v(reaction_i);
            end
        
            D = S*diag(v1)*S';
            J = S*dv';
            % solve J*C + C*J' + D = 0
            C = lyap(J,D);
        
            exp_CV(alter_i, level_i, reaction_i + 1)  = sqrt(C(4,4)/p.^2);
        end
    end
 end


save exp_CV.mat exp_CV
save exp_mean.mat exp_mean
