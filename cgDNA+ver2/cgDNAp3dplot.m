function [] = cgDNAp3dplot(Data,varargin)

%--------------------------------------------------------------------------
% cgDNA function: [] = cgDNA3dplot(myData)
%--------------------------------------------------------------------------
% This function plots the 3D reconstruction of the ground states as rigid
% bodies colored according to the sequence and the Crick Watson pairing.
%
% Input:
%
% MyData        structure array as output by the function cgDNA.m
%
% varargin      See Note 2 for possible options
%
%
% Output:
%
%   figure 1    plot of rigid bodies displayed as patches. Color scheme:
%               red A, blue T, green G, and yellow C. 
%
%   Note 1: If MyData contains M sequences, then M structures will be
%           displayed
%
%   Note 2: The optional arguments are : 
%       - RotDeg: possible values 0 or 1, default 1.
%         The rotation part of the coordinates can be displayed in
%         degrees or in rad/5. By default the unit are degrees
%
%       - Axis: possible values 'off' or 'on', default 'off'
%         If 'on' is chosen the xyz axis are displayed wiht grid ansd axis 
%         labels.       
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Submitted (2014). 
%
%--------------------------------------------------------------------------

% Check varargin for optional input
p = CheckInput(varargin);

% Add the folder containing the Viewer to the path
addpath(genpath('./cgDNApviewer'));
addpath('./Utilities');


% Initialise the figure
%figure 

% For loop over the number of sequences
for i = 1:length(Data)
    % Plot the 3d structrue of the ith sequence
    cgDNApviewer(Data(i).groundstate, Data(i).sequence,p.Reframe) ;
end

% By default do not plot axis
axis(p.flag_axis)

hold off

end

function p = CheckInput(inputarg)

p = inputParser ;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validUnit   = @(x) isnumeric(x) && x == 1 || x == 0;

addOptional(p,'Reframe',1,validScalarPosNum) ;
addOptional(p,'Axis',0,validUnit) ;


parse(p,inputarg{:});

p=p.Results;

if p.Axis == 0 
    p.flag_axis = 'off';
    
elseif p.Axis == 1
    p.flag_axis = 'on';
    
end

end

