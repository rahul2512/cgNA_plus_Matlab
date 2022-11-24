function Data = cgDNAp(varargin)

%--------------------------------------------------------------------------
% cgDNA function: myData = cgDNA(varargin) 
%--------------------------------------------------------------------------
%Constructs the cgDNA ground-state vector and stiffness matrix for one or a
%collection of DNA sequences.

% INPUT: can either take the form: 

%    1) cgDNA(seq, param_nbr, seq_label), with
%
%
%   seq             sequence along reference strand
%
%   param_nbr       integer (1,2,3 or 4) labelling the desired parameter set. 
%                   Default is 4.
%
%   seq_label       name (string) of the sequence (optional)

% OR 2) cgDNA(array), where array is a (n*k) cell array (k=1, 2 or 3), where each of the
%       n lines contains a sequence, a parameter set number (set to 4 by default
%       if k = 1), and a sequence label (optional).

%--------------------------------------------------------------------------

% OUTPUT: (1*n) array of Matlab struct MyData, with fields

%   seqlabel            name of the sequence
% 
%   sequence            sequence along reference strand
% 
%   nbp                 length of the sequence (number of basepairs)
%
%   paramset            name of the parameter set used
% 
% 
%   groundstate         ground-state coordinate vector 
%                       [size N x 1], 
%
%   stiff               stiffness matrix
%                       [size N x N],
%
%   where N = 12*nbp - 6 and nbp is the length 
%   of the sequence (number of basepairs). 


%--------------------------------------------------------------------------
% EXAMPLES


% 1) MyData = cgDNA('ATAGTA') OR MyData = cgDNA('ATAGTA', 4) returns a
%   struct MyData with fields:
% 
%        seqlabel: []
%        sequence: 'ATAGTA'
%             nbp: 6
%        paramset: 'cgDNAps4'
%     groundstate: [66×1 double]
%           stiff: [66×66 double]


% 2) MyPoly = cgDNA({'AAAA',2, 'polyA'; 'TTTTT',3, 'PolyT'});
% 
%   Then MyPoly(1) is a struct with fields:
% 
%        seqlabel: 'polyA'
%        sequence: 'AAAA'
%             nbp: 4
%        paramset: 'cgDNAps2'
%     groundstate: [42×1 double]
%           stiff: [42×42 double]

%--------------------------------------------------------------------------
%
%
%
%--------------------------------------------------------------------------

%If only one sequence as input, reformat input as a cell array
if nargin < 4 && ~iscell(varargin{1})
    for i=1:nargin
        input(i) = varargin(i);
    end
elseif nargin == 1 && iscell(varargin{1})
    input = varargin{1} ;
end
    
%Store sequence
sequence = input(:,1) ;

%Store parameter set number (default is cgDNA+ps1)
if size(input,2) < 2
    paramset_nbr = num2cell(ones(size(input,1),1));
else
    paramset_nbr = input(:,2) ;
end

%Store sequence labels (if any)
flag_seqlabel = 0 ;
if size(input,2) > 2
    seqlabel = input(:,3) ;
    flag_seqlabel = 1 ;
end
    

nbr_seq = length(sequence) ;

%Create output Data
Data = struct( ...
    'seqlabel' , cell(1,nbr_seq), ...
    'sequence', cell(1,nbr_seq) , ...
    'nbp', cell(1,nbr_seq) , ...
    'paramset', cell(1,nbr_seq) , ...
    'groundstate', cell(1,nbr_seq) ,  ...
    'stiff', cell(1,nbr_seq) ...
    );

%%Fill output for each sequence
for i = 1:nbr_seq
    if flag_seqlabel && isempty(seqlabel{i})~=1
        Data(i).seqlabel = seqlabel{i} ;
    else 
        Data(i).seqlabel = [ 'Seq' num2str(i) ] ;
    end
    
    Data(i).sequence = sequence{i} ;
    Data(i).nbp = length(sequence{i}) ; 
    Data(i).paramset = [ 'cgDNA+ps' num2str(paramset_nbr{i}) ] ;
    
    PS = load( [ 'ParameterSets/' Data(i).paramset ]) ;
    
    % Construct cgDNA shape vector and stiffness matrix
    [Data(i).groundstate, Data(i).stiff] = constructSeqParms(sequence{i}, PS) ;
    
end

end

