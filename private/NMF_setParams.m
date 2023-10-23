function [NMFparams,X_ft]=NMF_setParams(X_ft,param,NMFparams)
% [NMFparams,X_ft]=NMF_setParams(cfreq_amplitudes,param,NMFparams)
% Set NMF default parameters
%
% NMFparams.
%           NUM_MAX_ITER - Maximum number of NMF iterations (150)
%           ALPHA_A - Temporal smoothness value of A (50)
%           B - Beta divergence (0 itakura, 1 Kullback, 2 EUC) (1)
%           ETA_S - Update convergence speed for S (0.5)
%           ETA_A - Update convergence speed for A (0.5)
%           n_max - Excitation filter number (1)
%           m_max - number of partials (20)
% Input arguments: 
%	cfreq_amplitudes = spectral amplitude x time matrix
%   param = Sinusoidal model parameters
%   NMFparams* = Create if not exits-
%
% Output: 
%	NMFparams = NMF function parameters
%
% Julio Carabias / Francisco Rodriguez Diciembre 2011

% Initize struc
if nargin<3,
    NMFparams=struct();
end;

% NMF iterations
if ~isfield(NMFparams,'NUM_MAX_ITER')
    NMFparams.NUM_MAX_ITER=150;end;
% Gains temporal constraint
if ~isfield(NMFparams,'ALPHA_A')
    NMFparams.ALPHA_A=50;end;
% Sparsity Constraint
if ~isfield(NMFparams,'lambda')
    NMFparams.lambda=0;end;  
% Beta divergencia (0 itakura, 1 Kullback, 2 EUC)
if ~isfield(NMFparams,'B')
    NMFparams.B=1;end;
% control update S
if ~isfield(NMFparams,'ETA_S')
    NMFparams.ETA_S=0.5;end;
% control update A
if ~isfield(NMFparams,'ETA_A')
    NMFparams.ETA_A=0.5;end;  

% No instruments
if ~isfield(NMFparams,'j_max')
NMFparams.j_max=1;end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of harmonics
NMFparams.harmonics=unique(round(12*log2(1:NMFparams.m_max)));

% Format input matrix
if size(X_ft,1)>numel(param.midi_min:param.midi_max),
    X_ft = X_ft(param.midi_min:param.midi_max,:);
end;

% NMF params
NMFparams.f_max=size(X_ft,1);
NMFparams.t_max=size(X_ft,2); 
NMFparams.p_max=param.midi_max-param.midi_min+1;
NMFparams.p_min = param.midi_min;

return;