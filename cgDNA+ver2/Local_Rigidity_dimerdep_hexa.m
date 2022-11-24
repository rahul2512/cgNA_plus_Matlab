function res = Local_Rigidity_dimerdep_hexa()

addpath('./Utilities/')

load all_dimers_hexa_context.mat

pos = 7;

dimer ={'TT';'CC';'CT';'TC';
        'TA';'TG';'CG';'CG';
        'AT';'AC';'GC';'GT';
        'AA';'AG';'GA';'GG'};
 
    
for k = 1:16
    
    seqs = hexa.(dimer{k});
    
    for i = 1:256
        
        seq = [ 'GC' RndSeq(2) seqs(i,:) RndSeq(2) 'GC' ];
        Data = cgDNAp(seq);
        
        bp_level = frames(Data.groundstate);
        
        R = bp_level(pos).R;
        r = bp_level(pos).r; 
        
        Rup = bp_level(pos+1).R;
        rup = bp_level(pos+1).r;
        
        Rdown = bp_level(pos-1).R; 
        rdown = bp_level(pos-1).r;
        
        up(1:3,i) = diag(R'*Rup);
        down(1:3,i) = diag(R'*Rdown);
        
        up(4:6,i) = R'*(rup-r);
        down(4:6,i) = R'*(rdown-r);     
    end
    
    res.(dimer{k}).up = up ;
    res.(dimer{k}).down = down ;
    
%     subplot(4,4,k)
%     histogram(up(6,:),'BinWidth',0.05);
%     xlim([-4 4]);
%     ylim([0 40]);
%     hold on
%     histogram(down(6,:),'BinWidth',0.05);
%     title(dimer{k})
%     xlim([-4 4]);
%     ylim([0 40]);
    
end

save('Local_Rigidity.mat','res');

end


function seq = RndSeq(n)

base = {'A','T','C','G'};

seq = char(n);
id = randi(4,1,n);
for i = 1:n
    seq(i) = base{id(i)};
end

end

