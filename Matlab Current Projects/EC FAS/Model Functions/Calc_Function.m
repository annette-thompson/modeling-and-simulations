%% Calc_Function
function [FA_raw, rel_rate, NAD_rate] = Calc_Function(T,C,S)
% Calculates production profile and rate
%   Input:
%       T: time points [sec]
%       C: concentration matrix [uM]
%       S: structure storing all variables
%   Output:
%       FA_raw: final concentrations of all free fatty acids [uM]
%       rel_rate: rate of fatty acid production over total time [uM/min]
%       NAD_rate: conversion of NAD(P)H to NAD(P)+ over total time [uM/min]

FA_indices = contains(S.labels, '_FA', 'IgnoreCase', false);
FA_raw = C(end, FA_indices);
FA_weighted = sum(S.FA_dist / 16 .* FA_raw);
NAD_converted = sum(C(end, strcmp(S.labels, 'c_NADP') | strcmp(S.labels, 'c_NAD')));

% Calculate rates
rel_rate = FA_weighted / T(end) * 60;
NAD_rate = NAD_converted / T(end) * 60;
