function [minb,maxb] = fcn_baseline(lc,latency,len_base)

% fcn_baseline() - function which returns the start and the end of 
%                  baselines.
%
% Author: Fanny Grosselin, 2018
%                       
% Usage:
%   >> [minb,maxb] = fcn_baseline(lc,latency,len_base)
%
% Inputs:
%   lc       - [raw vector of integers] It contains sample points indicating
%              starts or ends of baselines. It is the output of the function
%              find_resp_marks().
%   latency  - ['pos'|'neg'] If it's 'pos', lc will be starts of baselines.
%              If it's 'neg', lc will be ends of baselines.
%   len_base - [char] Number written in char. It describes length of
%              baselines.
%
% Outputs :
%   minb - [raw vector of integers] It contains starts of baselines.
%   maxb - [raw vector of integers] It contains ends of baselines.
%
% See also:
%   pop_evokedinduced(), find_resp_marks().
%
 
EEG = evalin('base','EEG');

%% Start of baseline
%  -----------------
lc = lc;
 
 %% Latency of baseline
 %  -------------------
 if strcmp(latency,'pos')
     minb = lc; % start of baseline = lc
     maxb = [];
 else
     minb = [];
     maxb = lc; % end of baseline = lc
 end;
 
 %% Length of baseline
 %  ------------------
 factor_time = abs(abs(EEG.times(1)) - abs(EEG.times(2)));
 if isempty(maxb)
     maxb = minb + str2num(len_base)/factor_time;
 elseif isempty(minb)
     minb = maxb - str2num(len_base)/factor_time;
 end;


