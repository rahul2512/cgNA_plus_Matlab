load Deca_context_res.mat

dimer ={'AT';'TA';'CG';'GC';
        'AA';'TT';'CC';'GG';
        'CT';'TC';'GA';'AG';
        'GT';'TG';'CA';'AC'};
figure
for k = 1:16    
    
    
    subplot(4,4,k)
    histogram(minor{k},'BinWidth',0.1,'Normalization','pdf');
    xlim([8 20]);
    ylim([0 5]);
    hold on
    histogram(major{k},'BinWidth',0.1,'Normalization','pdf');
    title(dimer{k})
    xlim([8 20]);
    ylim([0 5]);
    
    
end