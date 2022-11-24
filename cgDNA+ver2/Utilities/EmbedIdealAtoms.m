function Data = EmbedIdealAtoms(Data)

for i = 1:length(Data)
    
    seq = Data(i).sequence;
   
    nbp = Data(i).nbp;
   
    bp_level = frames(Data(i).groundstate) ;
    
    [ AtomsBase, AtomsPho ] = Basepair2Atoms(bp_level, seq) ;
    
    Chain1_base = [AtomsBase.W];
    Chain1_pho  = [AtomsPho.W];
    
    Chain2_base = [AtomsBase.C]; Chain2_base = Chain2_base(nbp:-1:1) ;
    Chain2_pho  = [AtomsPho.C];  Chain2_pho = Chain2_pho(nbp-1:-1:1) ;
    
    Data.AtomsCoord = ArrangeAtoms(Chain1_base,Chain1_pho,Chain2_base,Chain2_pho,nbp) ;
    
end

end

function AtomsCoord = ArrangeAtoms(Chain1_base,Chain1_pho,Chain2_base,Chain2_pho,nbp)

% =========================================================================
%                               Chain 1
% =========================================================================

tmpCoord = zeros(3,1);
tmpNames = cell(1,1);

nbr_atoms = size(Chain1_base(1).coord,2);

tmpCoord(:,1) = Chain1_pho(1).coord(:,2);
tmpNames(1) = Chain1_pho(1).names(:,2);

tmpCoord(:,end+1:end+nbr_atoms) = Chain1_base(1).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain1_base(1).names(:,1:end);

AtomsCoord(1).W.coord = tmpCoord ;
AtomsCoord(1).W.names = tmpNames ;

for j=2:nbp-1
    
    tmpCoord = zeros(3,4);
    tmpNames = cell(1,4);
    
    nbr_atoms = size(Chain1_base(j).coord,2);
    
    tmpCoord(:,1:4) = Chain1_pho(j-1).coord(:,[1 4 5 3]);
    tmpNames(1:4) = Chain1_pho(j-1).names(:,[1 4 5 3]);
    
    tmpCoord(:,end+1) = Chain1_pho(j).coord(:,2);
    tmpNames(end+1) = Chain1_pho(j).names(:,2);
        
    tmpCoord(:,end+1:end+nbr_atoms) = Chain1_base(j).coord(:,1:end);
    tmpNames(end+1:end+nbr_atoms) = Chain1_base(j).names(:,1:end);
    
    AtomsCoord(j).W.coord = tmpCoord ;
    AtomsCoord(j).W.names = tmpNames ;
    
end

tmpCoord = zeros(3,4);
tmpNames = cell(1,4);

nbr_atoms = size(Chain1_base(nbp).coord,2);

tmpCoord(:,1:4) = Chain1_pho(nbp-1).coord(:,[1 4 5 3]);
tmpNames(1:4) = Chain1_pho(nbp-1).names(:,[1 4 5 3]);

tmpCoord(:,end+1:end+nbr_atoms) = Chain1_base(nbp).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain1_base(nbp).names(:,1:end);

AtomsCoord(nbp).W.coord = tmpCoord ;
AtomsCoord(nbp).W.names = tmpNames ;


% =========================================================================
%                               Chain 2
% =========================================================================

nbr_atoms = size(Chain2_base(1).coord,2);

tmpCoord = zeros(3,1);
tmpNames = cell(1,1);

tmpCoord(:,1) = Chain2_pho(1).coord(:,2);
tmpNames(1) = Chain2_pho(1).names(:,2);

tmpCoord(:,end+1:end+nbr_atoms) = Chain2_base(1).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain2_base(1).names(:,1:end);

AtomsCoord(1).C.coord = tmpCoord ;
AtomsCoord(1).C.names = tmpNames ;

for i=2:nbp-1
    
    tmpCoord = zeros(3,4);
    tmpNames = cell(1,4);
    
    nbr_atoms = size(Chain2_base(i).coord,2);
    
    tmpCoord(:,1:4) = Chain2_pho(i-1).coord(:,[1 4 5 3]);
    tmpNames(1:4) = Chain2_pho(i-1).names(:,[1 4 5 3]);
    
    tmpCoord(:,end+1) = Chain2_pho(i).coord(:,2);
    tmpNames(end+1) = Chain2_pho(i).names(:,2);
    
    tmpCoord(:,end+1:end+nbr_atoms) = Chain2_base(i).coord(:,1:end);
    tmpNames(end+1:end+nbr_atoms) = Chain2_base(i).names(:,1:end);
    
    AtomsCoord(i).C.coord = tmpCoord ;
    AtomsCoord(i).C.names = tmpNames ;
    
end

tmpCoord = zeros(3,4);
tmpNames = cell(1,4);

nbr_atoms = size(Chain2_base(nbp).coord,2);

tmpCoord(:,1:4) = Chain2_pho(nbp-1).coord(:,[1 4 5 3]);
tmpNames(1:4) = Chain2_pho(nbp-1).names(:,[1 4 5 3]);

tmpCoord(:,end+1:end+nbr_atoms) = Chain2_base(nbp).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain2_base(nbp).names(:,1:end);

AtomsCoord(nbp).C.coord = tmpCoord ;
AtomsCoord(nbp).C.names = tmpNames ;

end

