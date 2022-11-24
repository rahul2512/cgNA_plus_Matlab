function [shapes, stiff] = constructSeqParms(seq, paramset)
% -------------------------------------------------------------------------
% cgDNAp function: [shapes,stiff] = constructSeqParms(seq, paramset)
% -------------------------------------------------------------------------
% This function constructs the ground-state coordinate
% vector and stiffness matrix in non-dimensional Curves+
% form for a given sequence, using the specified parameter
% set in params.  
%
%
% Input: 
%
%   seq     sequence along reference strand;
%
%   params  a cgDNAplus parameter set. See Note 1 for details. 
%
%
% Output:
%
%   shapes  ground-state coordinate vector 
%           [size N x 1]
%
%   stiff   ground-state stiffness matrix
%           [size N x N]
%
%   where N = 12*nbp - 6 and nbp is the length 
%   of the sequence seq (number of basepairs). 
%
%
% Note 1: The cgDNAplus paramset is given as a MATLAB structure which
%         contains the following fields : 
%
%         sigma_int  : Contains the values of sigma for the iterior dimers
%         sigma_end5 : Contains the values of sigma for the 5' end dimers
%         sigma_end3 : Contains the values of sigma for the 3' end dimers
%
%         stiff_int  : Contains the stiff. blocks for the iterior dimers
%         stiff_end5 : Contains the stiff. blocks for the 5' end dimers
%         stiff_end3 : Contains the stiff. blocks for the 3' end dimers
%         
%         Each field contains itself one field corresponding to each dimer.
%         For the moment only GC ends blocks are available in the
%         paramseter set.
%
%
% Note 2: (TO BE MODIFIED WITH NAMES OF THE PHOSPHATE GROUP COORDINATES)
%
%    The entries in the variables shapes and stiff
%    are consistent with the following ordering of
%    the structural coordinates
%
%     y_1, z_1, ..., y_{nbp-1}, z_{nbp-1}, y_{nbp}
%
%    where for each a=1,2,3,... we have
%
%     y_a = (Buckle,Propeller,Opening,Shear,Stretch,Stagger)_a
%
%     z_a = (Tilt,Roll,Twist,Shift,Slide,Rise)_a.
%
%    For example
%
%     shapes((a-1)*12+1) = Buckle at basepair a
%     shapes((a-1)*12+2) = Propeller at basepair a
%      ...
%     shapes((a-1)*12+6) = Stagger at basepair a
%     shapes((a-1)*12+7) = Tilt at junction a, a+1
%     shapes((a-1)*12+8) = Roll at junction a, a+1
%      ...
%     shapes((a-1)*12+12) = Rise at junction a, a+1.
%
%    Correspondingly, we have
%
%     stiff(i,j) = stiffness coefficient for the pair of
%                  coordinates shapes(i) and shapes(j).
%
%
% If you find this code useful, please cite:
%
%
% -------------------------------------------------------------------------

% Initialize variables
seq = upper(seq);
nbp = numel(seq);
N =24*nbp-18 ;
stiff_cell = cell(nbp-1,1);
sigma = zeros(N,1);

dimer = repmat(seq, [nbp-1 1]) ; 
dimer = [diag(dimer),diag(dimer,1)] ;

%5' endblocks
[I,J,stiff_block] = find(paramset.stiff_end5.(dimer(1,:)) );
stiff_cell{1}= [I,J,stiff_block];
sigma(1:36) = paramset.sigma_end5.(dimer(1,:))  ;


%Stiffness Matrix and sigma vector's core
for i = 2:nbp-2
    k = (i-2)*24+18;
    
    %Sparsity pattern indices and corresponding dimer block values are
    %stored
    [I,J,stiff_block] =  find( paramset.stiff_int.(dimer(i,:)));
    stiff_cell{i}= [I+ k*ones(size(I)),J+ k*ones(size(J)),stiff_block];
    sigma(k+1:k+42,1) = sigma(k+1:k+42,1) + paramset.sigma_int.(dimer(i,:));
    
end

k = (nbp-3)*24 + 18;

%3' endblocks
[I,J,stiff_block] = find(paramset.stiff_end3.(dimer(nbp-1,:)));
stiff_cell{nbp-1}= [I + k*ones(size(I)),J + k*ones(size(J)),stiff_block];

sigma(k+1:k+36,1) = sigma(k+1:k+36,1)+ paramset.sigma_end3.(dimer(nbp-1,:)) ;


%Assembly of the Sparse Stiffness matrix
IJV= cell2mat(stiff_cell);
stiff = sparse(IJV(:,1),IJV(:,2),IJV(:,3),N,N);

shapes = stiff\sigma;

end