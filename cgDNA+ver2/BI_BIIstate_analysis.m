function res = BI_BIIstate_analysis(dimer)

addpath('./Utilities/')

load all_dimers_hexa_context

seqs = hexa.(dimer);

pos = 10;

for i = [1 6]
    
   
    seq = [ 'GC' RndSeq(4) seqs(i,:) RndSeq(4) 'GC' ];
    
    Data = cgDNAp(seq);
    Data = SugarRingAtoms(Data,pos) ;
    Data = BI_BIIstates(Data,pos) ;
    res{i} = Data.AtomsCoord(pos).W.state;
    %save( [ './Res_BI_BII_states/' dimer '.mat' ] , 'res' )

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
