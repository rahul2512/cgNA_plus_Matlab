function plot_Groove_RYalph()

dimer ={'AT';'TA';'CG';'GC';
    'AA';'TT';'CC';'GG';
    'CT';'TC';'GA';'AG';
    'GT';'TG';'CA';'AC'};

flank = {'YY','YR','RY','RR'};

load RY_alphabet_deca_data.mat %#ok<*LOAD>
load Deca_context_res.mat


res_mi = cell(4,4);
res_Mj = cell(4,4);

for k = 1:4
    
    res_tmp_mi = [];
    res_tmp_Mj = [];
    %for i = YY
     for i = 9   
        mi = minor{i};
        Mj = major{i};
        
        id = YY_f.(dimer{i}).(flank{k});
        res_tmp_mi = [res_tmp_mi mi(id)] ;
        res_tmp_Mj = [res_tmp_Mj Mj(id)] ;
        
    end
    
    res_mi{k,1} = res_tmp_mi;
    res_Mj{k,1} = res_tmp_Mj;
    
    res_tmp_mi = [];
    res_tmp_Mj = [];
   % for i = YR
     for i = 2     
        mi = minor{i};
        Mj = major{i};
        
        id = YR_f.(dimer{i}).(flank{k});
        res_tmp_mi = [res_tmp_mi mi(id)] ;
        res_tmp_Mj = [res_tmp_Mj Mj(id)] ;
        
    end
    
    res_mi{k,2} = res_tmp_mi;
    res_Mj{k,2} = res_tmp_Mj;
    
    res_tmp_mi = [];
    res_tmp_Mj = [];
  %  for i = RY
     for i = 4   
        mi = minor{i};
        Mj = major{i};
        
        id = RY_f.(dimer{i}).(flank{k});
        res_tmp_mi = [res_tmp_mi mi(id)] ;
        res_tmp_Mj = [res_tmp_Mj Mj(id)] ;
        
    end
    
    res_mi{k,3} = res_tmp_mi;
    res_Mj{k,3} = res_tmp_Mj;
    
    res_tmp_mi = [];
    res_tmp_Mj = [];
    
    
   % for i = RR
     for i = 12   
        mi = minor{i};
        Mj = major{i};
        
        id = RR_f.(dimer{i}).(flank{k});
        res_tmp_mi = [res_tmp_mi mi(id)] ;
        res_tmp_Mj = [res_tmp_Mj Mj(id)] ;
        
    end
    
    res_mi{k,4} = res_tmp_mi;
    res_Mj{k,4} = res_tmp_Mj;
    
end

figure
count = 1 ;
for i = 1:4
    
    for j = 1:4
        
        subplot(4,4,count)
        histogram(res_mi{i,j},'BinWidth',0.1,'Normalization','pdf');
        xlim([8 20]);
        ylim([0 5]);
        hold on
        histogram(res_Mj{i,j},'BinWidth',0.1,'Normalization','pdf');
        
        if i == 1
            title(flank{j})
        end
        
        xlim([8 20]);
        ylim([0 5]);
        
        count = count + 1 ;
    end
    
end


end

