function [] = cgDNApviewer( shape, sequence , bp_id )

%--------------------------------------------------------------------------
% cgDNA function: [] = cgDNAviewer( shape, sequence  )
%--------------------------------------------------------------------------
% 3D viewer of cgDNA ground-states, or in generale, 3D viewer of a coarse
% grained DNA fragment given in cgDNA-like internal coordinates. 
%
% Input : 
%
%  shapes    ground-state coordinate vector 
%            in non-dimensional Curves+ form (see Note 1)
%            [size N x 1].
%
%  sequence  DNA sequence in the {A,T,G,C} alphabet corresponding to the
%            cofiguration 'shapes'
%
%  bp_id     Index of the basepair to be used as center frame (i.e
%            baspair(bp_id).R = eye(3) and basepair(bp_id).r = zeros(3,1). 
%
% Note 1:
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
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Submitted (2014). 
%
%--------------------------------------------------------------------------

% Compute the absolute coordinates from the internal coordinates
bp_level = frames(shape);

% If the third argument is provided, change the reference frame
if nargin > 2
    if bp_id ~=1
        % Reframe the base and basepair frames with respect to the
        % reference basepair index
        bp_level = reframe(bp_level,bp_id);
    end
end

% Compute 
[ RB, RBc, Pw, Pc ] = sequence3DConf( bp_level, sequence ) ;

% Display rigid bodies
display3Dstructure(RB, RBc, Pw, Pc)

% Fix the point of view
view(3)

% Display the axis labels
xlabel('x')
ylabel('y')
zlabel('z')

% Visualize the grid 
grid on 

% Fix the axis properties 
axis equal on

end

