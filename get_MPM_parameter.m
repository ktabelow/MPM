function parameter = get_MPM_parameter(interpolation_kind)
% returns the number corresponding to the interpolation method to use in
% spm_slice_vol
%
% FORMAT X = get_MPM_parameter(interpolation_kind);
% INPUTS
% interpolation_kind   -  
% give the number of the interpolation method for the resampling:
%    nearest_neighbour                0  Zero-order hold (nearest neighbour).
%    trilinear_interpolation          1  First-order hold (trilinear interpolation).
%    lagrange                         2  Order Lagrange (polynomial) interpolation 
%    sinc                            -4  Order -4 of sinc interpolation.
% 
% CDA 02.06.2017
%

if nargin ==1
switch interpolation_kind
    case 'nearest_neighbour' 
        parameter = 0;
    case 'trilinear_interpolation'
        parameter = 1;
    case 'lagrange'
        parameter = 2;
    case 'sinc'
        parameter = -4; 
    otherwise
        parameter = 0;
end
else
    parameter = 0;
end