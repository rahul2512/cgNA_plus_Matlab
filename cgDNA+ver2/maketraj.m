function Data = maketraj(Data,filename,nbr_samples,FlagSugarAtoms)
%Rahul
%use maketraj(rand,'rand',200,0)
%Thus, given a data structre, this code generates a trajectory 
% of snapshots. The first snapshots is ground state and rest are Montecarlo
% samples.  
%This code first runs the monte_carlo_rahul thus save shape corresponinding
%to MC draws in config variables... Then use the slight modified version of
%makepdb (function is below not using the defualt alessandro's function)
% Also, not using the alessandro embedding function .. it is also modified
% and in the end of this file..
%Currently, not suitable for sugar fitting
Data = cgDNAp_MonteCarlo_rahul(Data,nbr_samples);
for i=0:nbr_samples
    makepdb(Data,filename,i,nbr_samples,FlagSugarAtoms);
end
end
%addpath(genpath('./Utilities'))
function Data = makepdb(Data,filename,index,nbr_samples,FlagSugarAtoms)

if nargin < 5
    FlagSugarAtoms = 0;
end

for l=1: length(Data)
    
    nbp = Data(l).nbp;
    
    if isfield(Data,'AtomsCoord')~=1
        
        if FlagSugarAtoms
            Data = SugarRingAtoms(Data) ;
        else
            Data = EmbedIdealAtoms(Data,index) ;
        end
    end
    % open sequence pdb file
    if length(Data) > 1 || (length(Data) == 1 && ~isempty(Data(l).seqlabel))
        if isempty(Data(l).seqlabel)
            fpdb = fopen([ 'PDBstruct/' filename '_Seq' num2str(l) '.pdb' ] ,'a+');
        else
            fpdb = fopen([ 'PDBstruct/' filename '_' Data(l).seqlabel '.pdb' ] ,'a+');
        end
    else
        fpdb = fopen(filename,'a+');
    end
    
    ntotal = 0;
    seq = Data(l).sequence;
    for i = 1:nbp %main strand
        
        Watson = [Data(l).AtomsCoord.W];
        
        if i == 1
            suffix = '5';
        elseif i == nbp
            suffix = '3';
        else
            suffix = [];
        end
        
        for j=1:length(Watson(i).names)
            ntotal = ntotal +1;
            fprintf(fpdb,'ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f\n', ntotal, Watson(i).names{j}, [ 'D' seq(i) suffix ], 'A', i, Watson(i).coord(:,j)');
        end
        
        
    end
    
    ntotal = ntotal +1;
    fprintf(fpdb,'TER\n');
    
    
    seq = wcc(Data(l).sequence,-1);
    for i = 1:nbp %complementary strand
        
        Crick = [Data(l).AtomsCoord.C];
        
        
        if i == 1
            suffix = '5';
        elseif i == nbp
            suffix = '3';
        else
            suffix = [];
        end
        
        
        
        for j=1:length(Crick(i).names)
            ntotal = ntotal +1;
            fprintf(fpdb,'ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f\n', ntotal, Crick(i).names{j}, [ 'D' seq(i) suffix ], 'B', i+nbp, Crick(i).coord(:,j)');
        end
        
    end
    
    ntotal = ntotal +1;
    if index < nbr_samples
        fprintf(fpdb,'TER\nENDMDL\n');
    else
        fprintf(fpdb,'TER\nENDMDL\nEND\n');
    end
    fclose(fpdb);
end

end


function Data = EmbedIdealAtoms(Data,index)

for i = 1:length(Data)
    seq = Data(i).sequence;
     
    nbp = Data(i).nbp;
    if index==0
        bp_level = frames(Data(i).groundstate) ;
    else
        bp_level = frames(Data(i).config(:,index)) ;
    end
    
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