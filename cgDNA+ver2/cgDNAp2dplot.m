function [] = cgDNAp2dplot(Data,varargin)

%--------------------------------------------------------------------------
% cgDNA function: [] = cgDNA2dplot(myData)
%--------------------------------------------------------------------------
% This function plots the ground-state coordinates
% to the screen.
%
% Input:
%
% MyData            structure array as output by the function cgDNA.m
%
% varargin          optional argument for the plot. See Note 3.
%
% Output:
%
%   panels 1...12   plot of each intra- and inter-
%                   basepair coordinate along the
%                   molecule
%
%
%
%   where N = 12*nbp - 6 and nbp is the length
%   of the sequence ( in number of basepairs).
%   Note that labeling for the graphics is optimised
%   for up to 20 or so base pairs.
%
% Note 1: If MyData contains M sequences, then M lines will
%         be drawn in each panels.
%
% Note 2: If MyData contains sequence label(s) they will be displayed
%
% Note 3: The optional arguments are :
%         - RotDeg: possible values 0 or 1, default 1.
%           The rotation part of the coordinates can be displayed in
%           degrees or in rad/5. By default the unit are degrees
%
%         - FontSize: possible values integers > 0, default 20.
%           Font size used for all the texts present in the plot (titles,
%           legend, axis labels, etc)
%
%         - LineWidth: possible values > 0, default 0.5.
%           Widht of all the plotted lines.
%
%        - Marker: possible values '.ox+*sdv^<>phn', default 'o'.
%           Marker to be used in the plot for all the sequences.
%
%         - MarkerkSize: possible values > 0, default 6.
%           Size of all the plotted marker. If myData contains less than 5
%           sequences no makers are used.
%
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks.
%  cgDNA: a software package for the prediction of sequence-dependent
%  coarse-grain free energies of B-form DNA. Submitted (2014).
%
%
%--------------------------------------------------------------------------
% Add Utilities to path
addpath('./Utilities');

% Check varargin for optional input
p = CheckInput(varargin);

% Read sequence label for the legend
lgtxt =  {Data.seqlabel};

% Find the longest sequences in the data structure
nbp = max([Data.nbp]) ;

% Initialise the 3d array that will contains the coordinates
nbr_seq = length(Data) ;
allCoordcgDNA = nan(nbp,12,nbr_seq) ;
allCoordPho = nan(nbp-1,12,nbr_seq) ;

% For loop over the number of sequence
for i  = 1:nbr_seq
    
    % Read ground state
    shape = Data(i).groundstate;
    
    % By default transfrom the groundstate in dimensional one
    if p.RotDeg
        shape = nondim2cur(shape);
    end
    
    % Decompose the shape into intra/inter and W/C pho rotations and
    % translations
    [Buckle_Propeller_Opening, ...
        Shear_Stretch_Stagger, ...
        pho_W_rot , pho_W_tr, ...
        Tilt_Roll_Twist, ...
        Shift_Slide_Rise, ...
        pho_C_rot , pho_C_tr ] = vector2shapes(shape);
    
    % Collects intras
    intra = [Buckle_Propeller_Opening, ...
        Shear_Stretch_Stagger];
    % Collects inters
    inter = [Tilt_Roll_Twist, ...
        Shift_Slide_Rise];
    % Collects Watson Phosphate
    pho_W = [pho_W_rot, ...
        pho_W_tr ] ;
    % Collects Crick Phosphate
    pho_C = [pho_C_rot, ...
        pho_C_tr ] ;
    
    % Store the decomposed cgDNA coordinates into 3 dimensional array
    tmp = [ [inter; NaN(1,6)] intra ];
    allCoordcgDNA(1:size(tmp,1),1:size(tmp,2),i) = tmp;
    
    % Store the decomposed Pho coordinates into 3 dimensional array
    tmp = [  [pho_W; NaN(1,6)] [ pho_C; NaN(1,6)] ];
    allCoordPho(1:size(tmp,1),1:size(tmp,2),i)=tmp;
    
end

% Initialise the figure
hsize = get(0,'ScreenSize');
fh=figure('Position',hsize);

% Define x axis limits and ticks
xlim = [0 nbp+1];
xtick = setdiff(sort([ 1 0:1:nbp]),0);
x = 1:nbp;

% For loop over the 12 panels for cgDNA coordinates
for j=1:12
    
    subplot(4,3,j);
    
    if j < 6
        dx = 0.5;
    elseif j > 6
        dx = 0;
    end
    
    % Read the data to plot in jth panel
    coor = squeeze(allCoordcgDNA(:,j,:));
    
    % Plot the data
    plot(x+dx, coor,'LineWidth', p.LineWidth,'Marker',p.Marker,'MarkerSize',p.MarkerSize);
    
    % Add title, x, and y label
    set(gca,'XLim',xlim,'XTick',xtick,'xticklabels',{'A','T','A','T','A','T','A','T','A','T','A','T/A','A','T','A','T','A','T',},...
        'Fontsize',p.FontSize);
    ylabel(p.ylabel{j}) ;
    xlabel(p.xlabel{j});
    title(p.title{j});
    
    grid on
end

% Show the legend
legend1=legend(lgtxt) ;
set(legend1,...
    'Position',[0.462211608886719 0.950662739322533 0.0901321411132813 0.0360824742268041],...
    'FontSize',p.FontSize);

% Set up axes position
axes('Position',[0,0,1,1],'visible','off');

hold off

% Initialise the figure
hsize = get(0,'ScreenSize');
fh=figure('Position',hsize);

% Define x axis limits and ticks
xlim = [0 nbp+1];
xtick = setdiff(sort([ 1 0:1:nbp]),0);
x = 1:nbp;

% For loop over the 12 panels for Pho coordinates
for j=1:12
    
    subplot(4,3,j);
    
    if j < 6
        dx = 1;
    elseif j > 6
        dx = 0;
    end
    
    % Read the data to plot in jth panel
    coor = squeeze(allCoordPho(:,j,:));
    
    % Plot the data
    plot(x+dx, coor,'LineWidth', p.LineWidth,'Marker',p.Marker,'MarkerSize',p.MarkerSize);
    
    % Add title, x, and y label
    set(gca,'XLim',xlim,'XTick',xtick,'xticklabels',{'A','T','A','T','A','T','A','T','A','T','A','T/A','A','T','A','T','A','T',},...
        'Fontsize',p.FontSize);
    ylabel(p.ylabel{j}) ;
    xlabel(p.xlabel_pho{j});
    title(p.title_pho{j});
    
    grid on
end

% Show the legend
legend1=legend(lgtxt) ;
set(legend1,...
    'Position',[0.462211608886719 0.950662739322533 0.0901321411132813 0.0360824742268041],...
    'FontSize',p.FontSize);

% Set up axes position
axes('Position',[0,0,1,1],'visible','off');

hold off

% Link all the panels for allowing simulatanous zooming
all_ha = findobj( fh, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' );

end

%--------------------------------------------------------------------------

function p = CheckInput(inputarg)

p = inputParser ;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validMarker = @(x) ischar(x) && contains('.ox+*sdv^<>phn',x);
validUnit   = @(x) isnumeric(x) && x == 1 || x == 0;

addOptional(p,'RotDeg',1,validUnit) ;
addOptional(p,'FontSize',20,validScalarPosNum) ;
addOptional(p,'LineWidth',0.5,validScalarPosNum) ;
addOptional(p,'Marker','o',validMarker) ;
addOptional(p,'MarkerSize',6,validMarker) ;

parse(p,inputarg{:});

p=p.Results;

if p.RotDeg
    
    p.ylabel = {
        'degrees'; 'degrees'; 'degrees'; ...
        'angstrom'; 'angstrom'; 'angstrom'; ...
        'degrees'; 'degrees'; 'degrees'; ...
        'angstrom'; 'angstrom'; 'angstrom' ...
        };
else
    p.ylabel = {
        'rad/5'; 'rad/5'; 'rad/5'; ...
        'angstrom'; 'angstrom'; 'angstrom' ; ...
        'rad/5'; 'rad/5'; 'rad/5'; ...
        'angstrom'; 'angstrom'; 'angstrom' ...
        };
end

p.title = {
    'Tilt'; 'Roll'; 'Twist'; ...
    'Shift'; 'Slide'; 'Rise'
    'Buckle'; 'Propeller'; 'Opening'; ...
    'Shear'; 'Stretch'; 'Stagger' ...
    };

p.title_pho = {
    'W Rot 1'; 'W Rot 2'; 'W Rot 3'; ...
    'W Tra 1'; 'W Tra 2'; 'W Tra 3' ; ...
    'C Rot 1'; 'C Rot 2'; 'C Rot 3'; ...
    'C Tra 1'; 'C Tra 2'; 'C Tra 3' ...
    };

p.xlabel = { 'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair'; 'basepair'; 'basepair'; ...
    'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair'; 'basepair'; 'basepair' ...
    };

p.xlabel_pho = { 'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair step'; 'basepair step'; 'basepair step'...
    };

end