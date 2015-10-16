% =========================================================================
% (c) Lars Ruthotto 2014/09/30
% http://www.mathcs.emory.edu/~lruthot/
%
% [I,omega,m,V] = loadImageSPM(PI,varargin)
%
% loads an nifti image from file
%
% Input:
%  PI     - filename
%
% Optional Arguments:
%  m      - resolution of target image, default: data resolution
%  Mref   - reference mat, default: data mat
%
% Output:
%  I      - image data, size(I) = m
%  omega  - domain size
%  m      - discretization size
%  V      - header of loaded image
%
% =========================================================================

function [I,omega,m,V] = loadImageSPM(PI,varargin)

if nargin==0,
    runMinimalExample;
    return;
end

m        = [];
Mref     = [];
res_hold = -7;
slices = [];

for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

V    = spm_vol(PI);
md   = V.dim;
if isempty(m), m = md; end;
if isempty(Mref), Mref = V.mat; end
if isempty(slices), slices = 1:m(3); end;
M    = V.mat;
tt   = sqrt(sum(M(1:3,1:3).^2)); 
omega = zeros(1,6);
omega(2:2:end) = tt(1:3).*md; 


I  = zeros([m(1:2), numel(slices)]);
M(:,1:3) = M(:,1:3)*diag(m./md);
Mt = Mref\ M;
for p=1:numel(slices)
    Ms = inv(spm_matrix([0 0 -slices(p) 0 0 0 1 1 1])*Mt);
    I(:,:,p) = spm_slice_vol(V,Ms,m(1:2),res_hold);
end
M(:,1:3) = M(:,1:3)*diag(md./m).^2;
m = size(I);
V.mat = M;
V.private.mat = M;
V.dim = m;


function runMinimalExample
m = [280,320,208];
PI = '/Home/stat/tabelow/projects/MPM/data/T1/s2013-10-31_14-54-151539-00001-00208-1.nii';
[I,omega,m1,V] = loadImageSPM(PI,'slices',[10,20,30]);
[I2,omega2,m2,V2] = loadImageSPM(PI,'slices', [10,20,30],'m',[ceil(m(1:2)./2) m(3)]);
[I3,omega3,m3,V3] = loadImageSPM(PI,'slices',[10,20,30],'m',[2*m(1:2) m(3)]);

figure(1); clf;
subplot(1,3,1)
imgmontage(I,omega,m1);

subplot(1,3,2)
imgmontage(I2,omega2,m2);
colormap gray

subplot(1,3,3)
imgmontage(I3,omega3,m3);
colormap gray
