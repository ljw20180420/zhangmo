function path = hmmMap(model, X)
% Find the most-probable (Viterbi) path through the HMM state trellis. 
%% Inputs:
% model - a struct as returned by e.g. hmmFit, which must contain
% at least the fields, pi and A for the starting state distribution
% and transition matrix respectively, as well as type, a string in {'gauss',
% 'discrete'}.
%
% X    - the local evidence vector (d-by-seqlen) where d=1 if type is
%        discrete.
%

% This file is from pmtk3.googlecode.com

pi = model.pi;
A  = model.A;
logB  = mkSoftEvidence(model.emission, X);
[path, j1, j2] = hmmViterbiC(log(pi+eps), log(A+eps), logB); %#ok<NASGU>
end
