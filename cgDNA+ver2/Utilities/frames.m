function bp_level = frames(shapes)  

%----------------------------------------------------------
% cgDNA function: [basepair] = frames(shapes)
%----------------------------------------------------------
% Given a ground-state coordinate vector in non-dimensional
% Curves+ form, this function constructs a reference point
% and frame for each base on each strand of the DNA
% according to the Tsukuba convention, and stores the
% result in the basepair structure (see Note 1). The
% reference point and frame vectors for each base are
% expressed relative to a fixed lab frame, which is assumed
% to coincide with the first basepair frame.
% 
%
% Input: 
%
%   shapes   ground-state coordinate vector 
%            in non-dimensional Curves+ form
%            [size N x 1].
%
%
% Output:
%
%   bp_levels    structure with reference point and frame
%                for each base and each phosphate group 
%                on each strand (see Note 1).
% 
%
% Note 1:
%
%   'bp_level' is a (1 x nbp) struct array with fields:
%    - 'D' : the frame of the base on the reading strand;
%    - 'r' : the coordinates of the base on the r. s.;
%    - 'Dc': the frame of the base on the complementary strand;
%    - 'rc': the coordinates of the base on the c. s.;
%
%    Reference point coordinates are 3x1 vectors, while frames 
%    are 3x3 matrices, with the frame coordinate vectors stored
%    as columns.  'nbp' is the length of the sequence.
%
%
% Note 2:
%
%    The entries of the input variable shapes must be 
%    consistent with the following ordering of the 
%    structural coordinates
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
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Submitted (2014). 
%
%----------------------------------------------------------

    nbp = (numel(shapes)+18)/24;

    % absolute coordinates of the first basepair 
    G = eye(3);
    q =[0, 0, 0]';
   
    % relative coordinates of the oligomer
    [eta, w, etapW, wpW, u, v, etapC, wpC] = vector2shapes(shapes);

    
    bp_level = InitializeStruct(nbp) ;
    
    for i=1:nbp
        
        bp_level(i).R = G ;
        bp_level(i).r = q ;
      
        % base pair:
        r = cay(eta(i,:));
        Gw = G * w(i,:)';
        
        % complimentary strand 
        bp_level(i).Rc = G * (sqrtm(r))'; 
        bp_level(i).rc = q - 0.5 * Gw;
        
        % main strand
        bp_level(i).Rw = bp_level(i).Rc * r;
        bp_level(i).rw = bp_level(i).rc + Gw;
        
        if i<nbp
            ru = cay(u(i,:));
            sqrtru = sqrtm(ru);
            H = G * sqrtru;
            % next base pair:
            G = G * ru;
            q = q + H * v(i,:)';
            
        end  

    end
    
    
    for i = 1:nbp - 1
    
      bp_level(i+1).Rpw = bp_level(i+1).Rw*cay(etapW(i,:)) ;
      bp_level(i+1).rpw  = bp_level(i+1).rw + bp_level(i+1).Rw*wpW(i,:)' ;
      
      Rc = bp_level(i).Rc*diag([1,-1,-1]) ;
      
      bp_level(i).Rpc = Rc*cay(etapC(i,:)) ;
      bp_level(i).rpc  = bp_level(i).rc + Rc*wpC(i,:)' ;
      
    end

end

%-------------------------------------------------------------

function bp_level = InitializeStruct(nbp)

bp_level  = struct('R', [],'r',[], ...
           'Rw', [], 'rw', [], ...
           'Rc', [], 'rc', [], ...
           'Rpw', [], 'rpw', [], ...
           'Rpc', [], 'rpc', cell(nbp,1)) ;


end

function Q = cay(k)

    I = eye(3) ;
    alpha = 1/10 ;
    k = alpha*k ;
    X = [   0   -k(3)  k(2) ;
           k(3)   0   -k(1) ;
          -k(2)  k(1)   0 ] ;
    Q = (I+X)/(I-X) ;

end
