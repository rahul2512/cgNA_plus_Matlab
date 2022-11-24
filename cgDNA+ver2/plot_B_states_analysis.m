function plot_B_states_analysis()

dimer ={'AT';'TA';'CG';'GC';
    'AA';'TT';'CC';'GG';
    'CT';'TC';'GA';'AG';
    'GT';'TG';'CA';'AC'};

flank = {'YY','YR','RY','RR'};

for k = 1:16
    
    load( [ 'Res_BI_BII_states/' dimer{k} '.mat' ]);
    
    Perc_I(k) = get_percentage_B_state(res);
    
end

Perc_II = 100-Perc_I;

figure
bar(1:16,[ Perc_I; Perc_II ]','stacked');
set(gca, 'XTickLabel',dimer, 'XTick',1:16)


load RY_alphabet_data.mat

res_matrix = nan(4,4);

row = 1;
for k = 1:4
    
    res_tmp = [];
    for i = YY
        
        load( [ 'Res_BI_BII_states/' dimer{i} '.mat' ]);      
        
        id = YY_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res(id)] ;
        
    end
    
    res_matrix(row,4) = get_percentage_B_state(res_tmp);
    
	row = row + 1; 
    
end

row = 1;
for k = 1:4
    
    res_tmp = [];
    for i = YR
        
        load( [ 'Res_BI_BII_states/' dimer{i} '.mat' ])      
        
        id = YR_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res(id)] ;
        
    end
    
    res_matrix(row,3) = get_percentage_B_state(res_tmp);
    
	row = row + 1; 
    
end

row = 1;
for k = 1:4
    
    res_tmp = [];
    for i = RY
        
        load( [ 'Res_BI_BII_states/' dimer{i} '.mat' ])
               
        
        id = RY_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res(id)] ;
        
    end
    
    res_matrix(row,2) = get_percentage_B_state(res_tmp);
    
    row = row + 1; 
    
end

row = 1;
for k = 1:4
    
    res_tmp = [];
    for i = RR
        
        load( [ 'Res_BI_BII_states/' dimer{i} '.mat' ])
        
        
        
        id = RR_f.(dimer{i}).(flank{k});
        res_tmp = [res_tmp res(id)] ;
        
    end
    
    res_matrix(row,1) = get_percentage_B_state(res_tmp);
    
    row = row + 1; 
    
end
figure
imagesc(res_matrix,[0 100])
title('BI state')

set(gca,'FontSize',50,'FontWeight','bold','Layer','top','XTickLabel',...
    {' ','RR',' ','RY',' ','YR',' ','YY',' '},'YTick',...
    [0.5 1 1.5 2 2.5 3 3.5 4 4.5],'YTickLabel',...
    {' ','Y..Y',' ','Y..R',' ','R..Y',' ','R..R',' '});

map=colormap('Bone');
map = map(end:-1:1,:);
colormap(map);


figure
imagesc(100-res_matrix,[0 100])
title('BII state')


set(gca,'FontSize',50,'FontWeight','bold','Layer','top','XTickLabel',...
    {' ','RR',' ','RY',' ','YR',' ','YY',' '},'YTick',...
    [0.5 1 1.5 2 2.5 3 3.5 4 4.5],'YTickLabel',...
    {' ','Y..Y',' ','Y..R',' ','R..Y',' ','R..R',' '});

map=colormap('Bone');
map = map(end:-1:1,:);
colormap(map);

end


function prec =  get_percentage_B_state(res)
count = 0 ;
for j = 1:length(res)
    
    switch res{j}
        
        case 'BI'
            count = count+1;
            
    end
    
end

prec = count / length(res) * 100 ;

end