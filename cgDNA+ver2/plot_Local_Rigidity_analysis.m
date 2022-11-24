function plot_Local_Rigidity_analysis(entry,xaxis,yaxis)

dimer ={'AT';'TA';'CG';'GC';
    'AA';'TT';'CC';'GG';
    'CT';'TC';'GA';'AG';
    'GT';'TG';'CA';'AC'};

flank = {'YY','YR','RY','RR'};

load Local_Rigidity.mat                                                         %#ok<*LOAD>


load RY_alphabet_data.mat

figure1 = figure ;

count = 1;
for k = 1:4
    
    res_tmp = [];                                                                    %#ok<*NASGU>
    for i = YY
                     
        id = YY_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res.(dimer{i}).up(:,id)] ;
        
    end

    RY_alph.YY.(flank{k}) = res_tmp;
    subplot(4,4,count)
    histogram(res_tmp(entry,:),'BinWidth',0.05)
    xlim(xaxis)
    ylim(yaxis)
    
    count = count + 1 ;
    
    res_tmp = [];
    for i = YR
                   
        id = YR_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res.(dimer{i}).up(:,id)] ;
        
    end
    
    RY_alph.YR.(flank{k}) = res_tmp;
    subplot(4,4,count)
    histogram(res_tmp(entry,:),'BinWidth',0.05)
    xlim(xaxis)
    ylim(yaxis)
   
	count = count + 1 ;
    
    res_tmp = [];
    for i = RY
                       
        id = RY_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res.(dimer{i}).up(:,id)] ;
        
    end
    
    RY_alph.RY.(flank{k}) = res_tmp;
    subplot(4,4,count)
    histogram(res_tmp(entry,:),'BinWidth',0.05)
    xlim(xaxis)
    ylim(yaxis)
    
	count = count + 1 ;
    
    res_tmp = [];
    for i = RR
                      
        id = RR_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res.(dimer{i}).up(:,id)] ;
        
    end
    
    subplot(4,4,count)
    histogram(res_tmp(entry,:),'BinWidth',0.05)
    xlim(xaxis)
    ylim(yaxis)
    
 	count = count + 1 ;
    
    RY_alph.RR.(flank{k}) = res_tmp;
    
end

annotation(figure1,'textbox',...
    [0.101 0.835481425322213 0.0157968749999999 0.0348749052312358],...
    'String','Y..Y',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.101 0.617892342683851 0.0157968749999999 0.0348749052312358],...
    'String','Y..R',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.10334375 0.392721758908264 0.0157968749999999 0.0348749052312358],...
    'String','R..Y',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.103734375 0.182714177407126 0.0157968749999999 0.0348749052312358],...
    'String','R..R',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.199828125 0.916603487490523 0.0157968749999999 0.0348749052312358],...
    'String','YY',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.40725 0.92191053828658 0.0157968749999998 0.0348749052312358],...
    'String','YR',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.615843750000001 0.91887793783169 0.0157968749999999 0.0348749052312358],...
    'String','RY',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.815062500000001 0.921152388172858 0.0157968749999998 0.0348749052312358],...
    'String','RR',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');



end