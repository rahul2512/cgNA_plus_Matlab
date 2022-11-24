function res = Local_Rigidity_dimerdep_hexa_2j()

addpath('./Utilities/')

load all_dimers_hexa_context.mat

pos = 7;

dimer ={'AT';'TA';'CG';'GC';
        'AA';'TT';'CC';'GG';
        'CT';'TC';'GA';'AG';
        'GT';'TG';'CA';'AC'};
 
    
for k = 1:16
    
    seqs = hexa.(dimer{k});
    
    for i = 1:256
        
        seq = [ 'GC' RndSeq(2) seqs(i,:) RndSeq(2) 'GC' ];
        Data = cgDNAp(seq);
        
        bp_level = frames(Data.groundstate);
        
        Rup = bp_level(pos+1).R;
        rup = bp_level(pos+1).r;
        
        Rdown = bp_level(pos-1).R; 
        rdown = bp_level(pos-1).r;
        
        res_tmp(1:3,i) = diag(Rdown'*Rup);       
        res_tmp(4:6,i) = Rdown'*(rup-rdown);
 
    end
    
    res.(dimer{k}).up = res_tmp ;
    
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

save('Local_Rigidity_2j.mat','res');

end


function seq = RndSeq(n)

base = {'A','T','C','G'};

seq = char(n);
id = randi(4,1,n);
for i = 1:n
    seq(i) = base{id(i)};
end

end

