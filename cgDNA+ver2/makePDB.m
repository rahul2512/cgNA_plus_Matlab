function Data = makePDB(Data,filename,FlagSugarAtoms)

%--------------------------------------------------------
% cgDNA function: [] = makePDB(Data, filename)
%--------------------------------------------------------
% Given the reference point and frame of each base,
% this function constructs the ideal coordinates of
% the non-hydrogen atoms of each base according to the
% Tsukuba definition, and writes the output to a PDB
% file (backbone atoms are not included).  The atomic
% coordinates are expressed relative to a fixed lab
% frame which coincides with the first basepair frame.
%
%
% Input:
%
%  Data             Data structure as output by cgDNAp.m
%
%  filename         Name of the PDB file to create
%
%  FlagSugarAtoms   If the flag is equal 1 compute the missing sugar rings
%                   atoms. If flag is zero just re-embed the ideal atoms
%                   for the base and the phosphate group.
%
% Auxiliary input:
%
%   idealBases.mat  matlab file with the ideal
%                   coordinates (in base frame) of
%                   the non-hydrogen atoms of the
%                   four bases T, A, C, G.
%
%
% Output: []
%
%
% Note 1:
%
%   'basepair' is a (1 x nbp) struct array with fields:
%    - 'Rw' : the frame of the base on the reading strand (Watson);
%    - 'rw' : the coordinates of the base on the r. s.;
%    - 'Rc': the frame of the base on the complementary strand (Crick);
%    - 'rc': the coordinates of the base on the c. s.;
%
%    Reference point coordinates are 3x1 vectors, while frames
%    are 3x3 matrices, with the frame coordinate vectors stored
%    as columns.  'nbp' is the length of the sequence.
%
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks.
%  cgDNA: a software package for the prediction of
%  sequence-dependent coarse-grain free energies of B-form
%  DNA. Nucleic Acids Research 2014; doi: 10.1093/nar/gku825.
%
%--------------------------------------------------------

addpath(genpath('./Utilities'))

if nargin < 3
    FlagSugarAtoms = 0;
end

for l=1: length(Data)
    
    nbp = Data(l).nbp;
    
    if isfield(Data,'AtomsCoord')~=1
        
        if FlagSugarAtoms
            Data = SugarRingAtoms(Data) ;
        else
            Data = EmbedIdealAtoms(Data) ;
        end
    end
    % open sequence pdb file
    if length(Data) > 1 || (length(Data) == 1 && ~isempty(Data(l).seqlabel))
        if isempty(Data(l).seqlabel)
            fpdb = fopen([ 'PDBstruct/' filename '_Seq' num2str(l) '.pdb' ] ,'w');
        else
            fpdb = fopen([ 'PDBstruct/' filename '_' Data(l).seqlabel '.pdb' ] ,'w');
        end
    else
        fpdb = fopen(filename,'w');
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
    fprintf(fpdb,'END');
    
    fclose(fpdb);
end

end
