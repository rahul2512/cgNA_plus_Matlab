function bp_level_new = reframe(bp_level,n)

%--------------------------------------------------------------------------
% cgDNA function: basepair_new = reframe(basepair,bp_id)
%--------------------------------------------------------------------------
% Given a structure basepair and an integer n, this function recompute a
% new basepair structure with basepair(n).R = eye(3) and basepair(n).r =
% zeros(3,1), i.e, this function center the absolute coordinates of the DNA
% fragment on the n basepair. If n = 1 the first basepair is the lab frame.
%
% Input:
%
%   basepair    structure with reference point and frame
%               for each base on each strand (see Note 1).
%
%   n           index of the basepair to use as center frame (see Note 2).
%               If not given, by default, the first basepair is used.
%
%
% Output:
%
%   basepair    structure with reference point and frame
%               for each base on each strand (see Note 1).
%
%
% Note 1:
%
%   'basepair' is a (1 x nbp) struct array with fields:
%    - 'R'  : the frame of the basepair;
%    - 'r'  : the coordinates of the basepair;
%    - 'Rw' : the frame of the base on the reading strand;
%    - 'rw' : the coordinates of the base on the r. s.;
%    - 'Rc' : the frame of the base on the complementary strand;
%    - 'rc' : the coordinates of the base on the c. s.;
%
%    Reference point coordinates are 3x1 vectors, while frames
%    are 3x3 matrices, with the frame coordinate vectors stored
%    as columns.  'nbp' is the length of the sequence.
%
% Note 2:
%
%   The center frame is actually the lab frames
%
%       O = eye(4)
%
%   If an index n is given (different from 1) all the frames
%   of the structure are moved according to a right rigid body
%   transformation given by the inverse of
%
%       g_n = |  R_n  r_n  |
%             | 0 0 0  1   |
%
%   where g_n is the nth basepair frames.
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks.
%  cgDNA: a software package for the prediction of sequence-dependent
%  coarse-grain free energies of B-form DNA. Submitted (2014).
%
%----------------------------------------------------------

nbp = length(bp_level);

if nbp < n
    errro('The basepair index given is invalid as is higher than the number of basepair');
end

% absolute coordinates of the first basepair
R_ref = bp_level(n).R;
r_ref = bp_level(n).r;

Q = R_ref';
q = -Q*r_ref;

bp_level_new = bp_level ;

P = diag([1,-1,-1]) ; 

for i=1:nbp
    
    bp_level_new(i).R = Q*bp_level(i).R;
    bp_level_new(i).r = Q*bp_level(i).r + q ;
    
    bp_level_new(i).Rw = Q*bp_level(i).Rw;
    bp_level_new(i).rw = Q*bp_level(i).rw + q ;
    
    bp_level_new(i).Rc = Q*bp_level(i).Rc;
    bp_level_new(i).rc = Q*bp_level(i).rc + q ;
    
    
    if i ~=1
        bp_level_new(i).Rpw = Q*bp_level(i).Rpw;
        bp_level_new(i).rpw = Q*bp_level(i).rpw + q ;
    end
    
    if i ~= nbp
        Rpc = bp_level(i).Rpc*P;
        
        bp_level_new(i).Rpc = (Q*Rpc)*P;
        bp_level_new(i).rpc = Q*bp_level(i).rpc + q ;
    end
    
end

end
