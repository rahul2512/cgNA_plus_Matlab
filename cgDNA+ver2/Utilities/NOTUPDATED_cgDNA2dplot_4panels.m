function [] = cgDNA2dplot_4panels(myData,varargin)

%--------------------------------------------------------------------------
% cgDNA function: [] = cgDNA2dplot_4panels(myData,varargin)
%--------------------------------------------------------------------------
% This function plots the ground-state coordinates to the screen divided
% into four panels : inter rotations, inter translation, intra rotations
% and intra translation
%
% Input:
%
% MyData            structure array as output by the function cgDNA.m
%
% varargin          optional argument for the plot. See Note 3.
%
% Output:
%
%   panels 1...4   plot of each intra- and inter-
%                   basepair coordinate along the
%                   molecule
%
%   where N = 12*nbp - 6 and nbp is the length
%   of the sequence ( in number of basepairs).
%   Note that labeling for the graphics is optimised
%   for up to 20 or so base pairs.
%
% Note 1: If MyData contains M sequences, then 3*M lines will
%         be drawn in each panels.
%
% Note 2: If MyData contains sequence label(s) they will be displayed with
%         the Line Style and Marker Style used in the plot.
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
%         - MarkerkSize: possible values > 0, default 6.
%           Size of all the plotted marker. If myData contains less than 5
%           sequences no makers are used. 
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Submitted (2014). 
%
%
%--------------------------------------------------------------------------

% Check varargin for optional input
p = CheckInput(varargin);

% Longest sequences in the data structure
nbp = max([myData.nbp]) ;

% Initialise the 3d array that will contains the coordinates
nbr_seq = length(myData) ;
allCoordcgDNA = nan(nbp,12,nbr_seq) ;

% Read the sequence labels
seqlegendtext = {myData.seqlabel};

% Initialise style variable 
LineStyle_id = 1;
MarkerStyle_id = 1;
PlotStyle = nan(nbr_seq,2);

for i  = 1:nbr_seq
    
    % Read ground state
    shape = myData(i).groundstate;
    
    % By default transfrom the groundstate in dimensional one
    if p.RotDeg
        shape = nondim2cur(shape);
    end
    
    % Decompose the shape into intra/inter rotations and translations
    [Buckle_Propeller_Opening, ...
        Shear_Stretch_Stagger, ...
        Tilt_Roll_Twist, ...
        Shift_Slide_Rise ] = vector2shapes(shape);
    
    % Collects intras
    intra = [Buckle_Propeller_Opening, ...
        Shear_Stretch_Stagger];
    % Collects inters
    inter = [Tilt_Roll_Twist, ...
        Shift_Slide_Rise];
    
    
    % Store the decomposed coordinates into 3 dimensional array
    tmp = [ [inter; NaN(1,6)] intra ] ;
    allCoordcgDNA(1:size(tmp,1),1:size(tmp,2),i) = tmp;
    
    % Store the Line and the Marker style
    PlotStyle(i,:) = [LineStyle_id , MarkerStyle_id] ;
    
    % Store the legned for the sequences
    seqlegendtext{i} = [ '  ' seqlegendtext{i} ' : ' ... 
                     ' Line ' p.LineStyle{LineStyle_id} ... 
                     ' Marker ' p.MarkerStyle{MarkerStyle_id}   '  '  ] ;
    
    % Only four line styles are available. For 5 or more sequence the line
    % style is combinied wiht a marker
    LineStyle_id = LineStyle_id + 1; 
    if LineStyle_id > 4 
    
        LineStyle_id = 1 ; 
        MarkerStyle_id  = MarkerStyle_id  + 1; 
        
    end
    
end

% Initialise the figure
hsize = get(0,'ScreenSize');
fh=figure('Position',hsize);

% Define x axis limits and ticks
xlim = [0 nbp+1];
xtick = setdiff(sort([ 1 0:5:nbp]),0);
x = 1:nbp;

% For loop over the 4 panels 
for j=1:4
    
    subplot(2,2,j);
    
    % Read the data to plot in jth panel
    range = 1+3*(j-1) : 3*j ;
    coor = squeeze(allCoordcgDNA(:,range,:));
    
    hold on 
    % Loop over all the sequences
    for i = 1:nbr_seq
        % Plot the data
        plot(x, squeeze(coor(:,:,i)),'LineStyle',p.LineStyle{PlotStyle(i,1)},'LineWidth', ...
            p.LineWidth,'Marker',p.MarkerStyle{PlotStyle(i,2)},'MarkerSize',p.MarkerSize);
        ax = gca;
        ax.ColorOrderIndex = 1;
                 
    end
    hold off
    
    % Add title, x, and y label 
    set(ax,'XLim',xlim,'XTick',xtick,...
        'Fontsize',p.FontSize);
    ylabel(p.ylabel{j}) ;
    title(p.title{j});
    
    % Show the legend
    legend1=legend(p.legend{j,:}) ;
    set(legend1,'FontSize',p.FontSize)
    
    grid on
    
end

hold off

% Write a super title containing the legend for the different sequences
a = axes;
t1 = title(cell2mat(seqlegendtext),'FontSize',p.FontSize,'FontWeight','normal');
a.Visible = 'off';
t1.Visible = 'on'; 

% Link all the panels for allowing simulatanous zooming
all_ha = findobj( fh, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' );


end

%--------------------------------------------------------------------------

function p = CheckInput(inputarg)

p = inputParser ;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validUnit   = @(x) isnumeric(x) && x == 1 || x == 0;

addOptional(p,'RotDeg',1,validUnit) ;
addOptional(p,'FontSize',20,validScalarPosNum) ;
addOptional(p,'LineWidth',0.5,validScalarPosNum) ;
addOptional(p,'MarkerSize',6,validScalarPosNum) ;


parse(p,inputarg{:});

p=p.Results;

if p.RotDeg
    
    p.ylabel = {
        'degrees'; 'angstrom' ;
        'degrees'; 'angstrom'
        };
else
    p.ylabel = {
        'rad/5'; 'angstrom';
        'rad/5'; 'angstrom'
        };
end

p.title = {
    'Inter Rotations'; 'Inter Translations'
    'Intra Rotations'; 'Intra Translations'
    };

p.legend = {'Tilt','Roll','Twist' ; 'Shift' 'Slide','Rise' ; ...
    'Buckle','Propeller','Opening'; 'Shear', 'Stretch' , 'Stagger'  };

p.LineStyle = {'-','--','-.',':'} ;
p.MarkerStyle = {'none';'o';'x';'+';'*';'s';'d';'v';'^';'<';'>';'p';'h';'n';'.'} ;

end