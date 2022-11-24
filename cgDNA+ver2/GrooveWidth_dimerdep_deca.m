function GrooveWidth_dimerdep_deca()

addpath('./Utilities/')

pos = 12;
dist =5;

dimer ={'AT';'TA';'CG';'GC';
        'AA';'TT';'CC';'GG';
        'CT';'TC';'GA';'AG';
        'GT';'TG';'CA';'AC'};

D=load('all_dimers_deca_context.mat');
 
parfor k = 1:16
   
    [minor{k} , major{k} ] = compute_groove(dimer{k},D,pos,dist);
    
end

save('Deca_context_res.mat', 'minor','major')

end

function [minor , major ] = compute_groove(dimer,D,pos,dist)

    res = nan(2*dist,256);
    seqs = D.deca.(dimer);
   
    for i = 1:length(seqs)
        
        seq = [ 'GC' RndSeq(4) seqs(i,:) RndSeq(4) 'GC' ];
        Data = cgDNAp(seq);
        
        bp_level = frames(Data.groundstate);
        
        Pw = bp_level(pos).rpw;
        
        for j = 1:dist
            
            Pc_up   = bp_level(pos+j-1).rpc;
            Pc_down = bp_level(pos-j).rpc;
            
            res(dist+j,i) = norm(Pc_up-Pw);
            res(dist-j+1,i) = norm(Pc_down-Pw);
            
        end
        
        minor(i) = min(res(dist:end,i));
        major(i) = min(res(1:dist-1,i));
        
    end

end

function seq = RndSeq(n)

base = {'A','T','C','G'};

seq = char(n);
id = randi(4,1,n);
for i = 1:n
    seq(i) = base{id(i)};
end

end

