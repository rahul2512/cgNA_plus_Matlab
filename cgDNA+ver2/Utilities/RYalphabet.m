function  RYalphabet()

load all_dimers_hexa_context

dimer ={'AT';'TA';'CG';'GC';
    'AA';'TT';'CC';'GG';
    'CT';'TC';'GA';'AG';
    'GT';'TG';'CA';'AC'};

YY = [ 7 9 10 6  ] ;
YR = [ 15 2 14 3 ] ;
RY = [ 4 13 1 16 ] ;
RR = [ 8 11 12 5 ] ;


YY_f_tmp.RR = [];
YY_f_tmp.RY = [];
YY_f_tmp.YR = [];
YY_f_tmp.YY = [];

RY_f_tmp.RR = [];
RY_f_tmp.RY = [];
RY_f_tmp.YR = [];
RY_f_tmp.YY = [];

YR_f_tmp.RR = [];
YR_f_tmp.RY = [];
YR_f_tmp.YR = [];
YR_f_tmp.YY = [];

RR_f_tmp.RR = [];
RR_f_tmp.RY = [];
RR_f_tmp.YR = [];
RR_f_tmp.YY = [];

flank = {'YY','YR','RY','RR'};

for k = YY
    
    seqs=hexa.(dimer{k});
    
    for i = 1:256
        
        h =seqs(i,:);
        
        switch h(2)
            
            case {'C','T'}
                type = 'Y';
                
            case {'G','A'}
                type = 'R';
        end
        
        switch h(5)
            
            case {'C','T'}
                type = [ type 'Y' ];
                
            case {'G','A'}
                type = [ type 'R' ];
        end
        
        YY_f_tmp.(type) = [ YY_f_tmp.(type) i ] ;
        
    end
       
    
    for i = 1:4
        YY_f.(dimer{k}).(flank{i}) = YY_f_tmp.(flank{i}) ;
    end
    
    
end

for k = YR
    
    seqs=hexa.(dimer{k});
    
    for i = 1:256
        
        h =seqs(i,:);
        
        switch h(2)
            
            case {'C','T'}
                type = 'Y';
                
            case {'G','A'}
                type = 'R';
        end
        
        switch h(5)
            
            case {'C','T'}
                type = [ type 'Y' ];
                
            case {'G','A'}
                type = [ type 'R' ];
        end
        
        YR_f_tmp.(type) = [ YR_f_tmp.(type) i ] ;
        
    end
        
    for i = 1:4
        YR_f.(dimer{k}).(flank{i}) = YR_f_tmp.(flank{i}) ;
    end
    
end

for k = RY
    
    seqs=hexa.(dimer{k});
    
    for i = 1:256
        
        h =seqs(i,:);
        
        switch h(2)
            
            case {'C','T'}
                type = 'Y';
                
            case {'G','A'}
                type = 'R';
        end
        
        switch h(5)
            
            case {'C','T'}
                type = [ type 'Y' ];
                
            case {'G','A'}
                type = [ type 'R' ];
        end
        
        RY_f_tmp.(type) = [ RY_f_tmp.(type) i ] ;
        
    end
        
    for i = 1:4
        RY_f.(dimer{k}).(flank{i}) = RY_f_tmp.(flank{i}) ;
    end
    
    
end

for k = RR
    
    seqs=hexa.(dimer{k});
    
    for i = 1:256
        
        h =seqs(i,:);
        
        switch h(2)
            
            case {'C','T'}
                type = 'Y';
                
            case {'G','A'}
                type = 'R';
        end
        
        switch h(5)
            
            case {'C','T'}
                type = [ type 'Y' ];
                
            case {'G','A'}
                type = [ type 'R' ];
        end
        
        RR_f_tmp.(type) = [ RR_f_tmp.(type) i ] ;
        
    end
    
    for i = 1:4
        RR_f.(dimer{k}).(flank{i}) = RR_f_tmp.(flank{i}) ;
    end
    
end

save('RY_alphabet_data.mat','YY','YR','RY','RR','YY_f','YR_f','RY_f','RR_f')

end

