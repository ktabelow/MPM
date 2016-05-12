% calculate P2_a and P2_b
%
% =========================================================================
% 2016/05/12
%
% adapted for MPM by C. D'Alonzo 
% original code in toolbox QVB
%
% 
%
% [P2_a,P2_b] = getPolynomsP2_ab(TR_pdw,TR_t1w,fa_pdw, fa_t1w)
% 
% calculate the coefficients of the 2nd degree polynoms P2_a and P2_b to be
% used in the correction of R1.
% Given a b1map correction field, these two polynoms have to be used so:
% A = P2_a(1)*b1map.^2+P2_a(2)*b1map+P2_a(3)
% B = P2_b(1)*b1map.^2+P2_b(2)*b1map+P2_b(3)
% and applyed as follow to calculate R1_corrected given R1_original
% R1_corrected= R1_original./(A*R1_original+B) (see calculateQI.m for this)
% 
%
% Input:
%
%  TR_pdw     - repetition time for PD weighted volumes
%  TR_t1w     - repetition time for T1 weighted volumes
%  fa_pdw     - flip angle for PD weighted volumes
%  fa_t1w     - flip angle for T1 weighted volumes
% 
%
% Output:
%
%
%  P2_a        - coefficients for P2_a polynom - 3 element-array
%  P2_b        - coefficients for P2_b polynom - 3 element-array
%
% =========================================================================


function [P2_a,P2_b]= getPolynomsP2_ab(TR_pdw,TR_t1w,fa_pdw, fa_t1w)

if nargin<4,
    error('There have to be 4 inputs: pdTR,t1TR,pdFA,t1FA');    
end

% Settings for R.Deichmann steady state correction using T2=64ms at 3T
% Correction parameters were calculated for 3 different parameter sets:
if (feq(TR_pdw, 23.7) && feq(TR_t1w, 18.7) && feq(fa_pdw, 6) && feq(fa_t1w, 20))
    % 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
    % PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
    disp('Classic FIL protocol');
    P2_a= [78.9228195006542,-101.113338489192,47.8783287525126];
    P2_b=[-0.147476233142129,0.126487385091045,0.956824374979504];    
elseif (feq(TR_pdw, 24.5) && feq(TR_t1w, 24.5) && feq(fa_pdw, 5) && feq(fa_t1w, 29))
    % 2) new FIL/Helms protocol
    % PD-weighted: TR=24.5ms; a=5deg; T1-weighted: TR=24.5ms; a=29deg
    disp('New FIL/Helms protocol');
    P2_a= [93.455034845930480,-120.5752858196904,55.911077913369060];
    P2_b=[-0.167301931434861,0.113507432776106,0.961765216743606];
elseif (feq(TR_pdw, 24.0) && feq(TR_t1w, 19.0) && feq(fa_pdw, 6) && feq(fa_t1w, 20))
    % 3) Siemens product sequence protocol used in Lausanne (G Krueger)
    %PD-weighted: TR=24ms; a=6deg; T1-weighted: TR=19ms; a=20deg
    disp('Siemens product Lausanne (GK) protocol');
    P2_a= [67.023102027100880,-86.834117103841540,43.815818592349870];
    P2_b=[-0.130876849571103,0.117721807209409,0.959180058389875];
elseif (feq(TR_pdw, 23.7) && feq(TR_t1w, 23.7) && feq(fa_pdw, 6) && feq(fa_t1w, 28))
    % 4) High-res (0.8mm) FIL protocol:
    % PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=23.7ms; a=28deg
    disp('High-res FIL protocol');
    P2_a= [1.317257319014170e+02,-1.699833074433892e+02,73.372595677371650];
    P2_b=[-0.218804328507184,0.178745853134922,0.939514554747592];
elseif (feq(TR_pdw, 25.25) && feq(TR_t1w, 25.25) && feq(fa_pdw, 5) && feq(fa_t1w, 29))
    % 4)NEW  High-res (0.8mm) FIL protocol:
    % PD-weighted: TR=25.25ms; a=5deg; T1-weighted: TR=TR=25.25ms; a=29deg
    disp('High-res FIL protocol');
    P2_a= [88.8623036106612,-114.526218941363,53.8168602253166];
    P2_b=[-0.132904017579521,0.113959390779008,0.960799295622202];
elseif (feq(TR_pdw, 24.5) && feq(TR_t1w, 24.5) && feq(fa_pdw, 6) && feq(fa_t1w, 21))
    % 5)NEW  1mm protocol - seq version v2k:
    % PD-weighted: TR=24.5ms; a=6deg; T1-weighted: TR=24.5ms; a=21deg
    disp('v2k protocol');
    P2_a= [71.2817617982844,-92.2992876164017,45.8278193851731];
    P2_b=[-0.137859046784839,0.122423212397157,0.957642744668469];
elseif (feq(TR_pdw, 25.0) && feq(TR_t1w, 25.0) && feq(fa_pdw, 6) && feq(fa_t1w, 21))
    % 6) 800um protocol - seq version v3* released used by MEG group:
    % TR = 25ms for all volumes; flipAngles = [6, 21 deg] for PDw and T1w
    % Correction parameters below were determined via Bloch-Torrey 
    % simulations but end result agrees well with EPG-derived correction 
    % for this RF spoiling increment of 137 degrees.
    % See: Callaghan et al. ISMRM, 2015, #1694
    disp('v3* 0.8mm R4 protocol');
    P2_a = [57.427573706259864, -79.300742898810441,  39.218584751863879];
    P2_b = [-0.121114060111119,   0.121684347499374,   0.955987357483519];
else
    disp('Spoiling correction not defined for this protocol. No correction being applied.');
    P2_a = [0, 0, 0];
    P2_b = [0, 0, 1];
end


function bl = feq(val, comp_val)
% floating point comparison
bl = abs(val - comp_val) <= eps(comp_val);
return;