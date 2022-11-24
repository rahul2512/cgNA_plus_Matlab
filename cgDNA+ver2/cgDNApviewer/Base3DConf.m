function [RB] = Base3DConf(r,R, base)

%--------------------------------------------------------------------------
% cgDNA function: [RB] = Base3DConf(r,R, base)
%--------------------------------------------------------------------------
% This function contructs a parallelepiped which rotation and position
% correspond to the rigid body configuration
%
%         RB = |    R  r |
%              | 0 0 0 1 |
% 
% and the color respect the base type given as argument.
%
% Input:
%
% r             rigid body's position.
%
% R             rigid body's orientation.
%
% base          base type corresponding to the rigid body configuration [R,r].
%               The base should be one of the standart alphabet {A,T,G,C}.
%               In order to construct rigid bodies for bases on the
%               complementary strand an extra 'c' should be added to base,
%               e.g., {'Ac','Tc','Gc','Cc'}.
%
%
% Note 1:The parallelepiped is colored using the color convention for the
% DNA bases: red A, blue T, green C, and yellow C
%
% Output:
%
% RB            structure containing two fields:
%               conf        verteces of the rigid body [ 8x3 matrix ]
%               color       the color of the rigid body [ 1x3 vector ]
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Submitted (2014). 
%
%--------------------------------------------------------------------------

% Define the dimension of the parallelepiped in the x and z direction
dz = 0.36 ;
dx = 4.2 ;
% Define the off set from the base frame
Y  = - 1.2 ;

% Depending on the base type, assign color and dimension on the y direction
switch base(1)
    case 'A'
        RB.color = [1 0 0] ;
        dy = 6 ; 
    case 'T'
        RB.color = [0 0 1] ;
        dy = 4.2 ;
    case 'G'
        RB.color = [0 1 0] ;
        dy = 6 ;
    case 'C'
        RB.color = [1 1 0] ;
        dy = 4.2 ;
end

% If complementary strand invert the sign on the y axis
if length(base) == 2
    dy = -dy ; 
    Y = -Y ;
end

% Define the vertices of the paralellepiped according to the base type
% around the origin
verteces_RB = [
    [0 0 0]    - 0.5*[ dx Y dz ]
    [dx 0 0]   - 0.5*[ dx Y dz ]
    [dx dy 0]  - 0.5*[ dx Y dz ]
    [0 dy 0]   - 0.5*[ dx Y dz ]
    [0 0 dz]   - 0.5*[ dx Y dz ]
    [dx 0 dz]  - 0.5*[ dx Y dz ]
    [dx dy dz] - 0.5*[ dx Y dz ]
    [0 dy dz]  - 0.5*[ dx Y dz ]
    ];

% Rotate and translate the parallelepiped according to base frame
RB.conf = (R*verteces_RB')' + repmat(r',[8,1]);

end

