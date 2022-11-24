dimer ={'AT';'TA';'CG';'GC';'TT';'CC';'GG';'CT';'TC';'GA';'AG';'GT';'TG';'CA';'AC'};

addpath('./Utilities/')

parfor k = 1:15

    BI_BIIstate_analysis(dimer{k});

    fprintf('dimer %s done \n', dimer{k})
    
end