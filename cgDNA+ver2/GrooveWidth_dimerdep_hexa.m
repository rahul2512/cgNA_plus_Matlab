function GrooveWidth_dimerdep_hexa()

addpath('./Utilities/')

load all_dimers_hexa_context.mat

pos = 12;
dist =5;

res = nan(2*dist,256);

dimer ={'AT';'TA';'CG';'GC';
        'AA';'TT';'CC';'GG';
        'CT';'TC';'GA';'AG';
        'GT';'TG';'CA';'AC'};

for k = 1:16
    
    seqs = hexa.(dimer{k});
    
    for i = 1:256
        
        seq = [ 'GC' RndSeq(6) seqs(i,:) RndSeq(6) 'GC' ];
        Data = cgDNAp(seq);
        
        bp_level = frames(Data.groundstate);
        
        Pw = bp_level(pos).rpw;
        
        for j = 1:dist
            
            Pc_up   = bp_level(pos+j-1).rpc;
            Pc_down = bp_level(pos-j).rpc;
            
            res(dist+j,i) = norm(Pc_up-Pw);
            res(dist-j+1,i) = norm(Pc_down-Pw);
            
        end
        
        minor(i,k) = min(res(dist:end,i));
        major(i,k) = min(res(1:dist-1,i));
        
    end
    
    subplot(4,4,k)
    histogram(minor,'BinWidth',0.1);
    xlim([9.5 18.5]);
    ylim([0 120]);
    hold on
    histogram(major,'BinWidth',0.1);
    title(dimer{k})
    xlim([9.5 18.5]);
    ylim([0 120]);
    
    save('Groove_width_hexa_res.mat','minor','major')
    
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

