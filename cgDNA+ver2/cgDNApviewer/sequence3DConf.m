function [RB,RBc,Pw,Pc] = sequence3DConf( bp_level, sequence )

%--------------------------------------------------------------------------
% cgDNA function: [RB , RBc] = sequence3DConf( basepair, sequence ) 
%--------------------------------------------------------------------------
% This function constructes two chains of parallelepipeds corresponding to
% the reading and compelmentary bases stored in basepair. The color of the
% boxes correspond to the given sequence.
%
% Input:
%
% bp_level   structure with reference point and frame
%                   for each base on each strand (see Note 1).
% 
% sequence   sequence along reference strand
%
%
% Output:
% 
% RB         structure containing two fields (reading strand):
%               conf        verteces of the rigid body [ 8x3 matrix ]
%               color       the color of the rigid body [ 1x3 vector ]
%
% RBc        structure containing two fields (complementary strand):
%               conf        verteces of the rigid body [ 8x3 matrix ]
%               color       the color of the rigid body [ 1x3 vector ]
%
%
% Note 1:
%
%   'basepair' is a (1 x nbp) struct array with fields:
%    - 'R' : the frame of the basepair;
%    - 'r' : the coordinates of the basepair;
%    - 'Rw' : the frame of the base on the reading strand;
%    - 'rw' : the coordinates of the base on the r. s.;
%    - 'Rc': the frame of the base on the complementary strand;
%    - 'rc': the coordinates of the base on the c. s.;
%
%    Reference point coordinates are 3x1 vectors, while frames 
%    are 3x3 matrices, with the frame coordinate vectors stored
%    as columns.  'nbp' is the length of the sequence.
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Submitted (2014). 
%
%--------------------------------------------------------------------------

% Store the total number of base pairs
nbp = length(sequence);

% Initialize the structures RB and RBc
RB = struct('color',cell(1,nbp),'conf',cell(1,nbp)) ;
RBc = struct('color',cell(1,nbp),'conf',cell(1,nbp)) ;

% For loop over the base pairs
for i = 1 : nbp
  
  % store the positions of the reading and complementary bases
  rw = bp_level(i).rw ;
  rc = bp_level(i).rc ;
  
  % store the orientation of the reading and complementary bases
  Rw = bp_level(i).Rw ;
  Rc = bp_level(i).Rc ;
  
  base = sequence(i) ;
  
  % Find the complementary base 
  switch base
    case 'A'
      baseC = 'Tc' ;
    case 'T'
      baseC = 'Ac' ;
    case 'G'
      baseC = 'Cc' ;
    case 'C'
      baseC = 'Gc' ;
  end
  
  % Compute the parallelepided corresponding to the base [Rw,rw]
  RB(i)  = Base3DConf(rw, Rw, base) ;
  
  % Compute the parallelepided corresponding to the base [Rw,rw]
  RBc(i) = Base3DConf(rc, Rc, baseC) ;
  
end

scale = 1.5 ;

nbp = length(bp_level) ;

Pc = zeros(4,3,nbp-1);
Pw = zeros(4,3,nbp-1);

for i = 1:nbp-1
    
    D = bp_level(i+1).Rpw ;
    O = bp_level(i+1).rpw ;
    Pw(:,:,i) = [O' ; (O+scale*D(:,1))' ; (O+scale*D(:,2))' ; (O+scale*D(:,3))'] ;
       
    D = bp_level(i).Rpc ;
    O = bp_level(i).rpc ;
    Pc(:,:,i) = [O' ; (O+scale*D(:,1))' ; (O+scale*D(:,2))' ; (O+scale*D(:,3))'] ;
     
end


end

