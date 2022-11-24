function [] = display3Dstructure(RB,RBc,Pw,Pc)

%--------------------------------------------------------------------------
% cgDNA function: [] = display3Dstructure(RB,RBc)
%--------------------------------------------------------------------------
% This script plots the DNA structure from the structures RB and RBc.
%
% Input:
%
% RB         structure containing two fields (reading strand):
%               conf        verteces of the rigid body [ 8x3 matrix ]
%               color       the color of the rigid body [ 1x3 vector ]
%
%
% RBc        structure containing two fields (complementary strand):
%               conf        verteces of the rigid body [ 8x3 matrix ]
%               color       the color of the rigid body [ 1x3 vector ]
%
%
% Output:
% figure 1   3D plot of the parallelepipeds stored in RB and RBc
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks.
%  cgDNA: a software package for the prediction of sequence-dependent
%  coarse-grain free energies of B-form DNA. Submitted (2014).
%
%--------------------------------------------------------------------------

% Define the faces of the rigid body
faces_matrix = [1 2 6 5
    2 3 7 6
    3 4 8 7
    4 1 5 8
    1 2 3 4
    5 6 7 8];

% Compute the total number of basepairs
nbp = length(RB);

% Plot the chains corresponding to the bases on the reading strand (RB)
% and on the complementary strand (RBC)
for j = 1 : nbp
    
    patch( 'Vertices',RB(j).conf, ...
        'Faces',faces_matrix,...
        'FaceVertexCData',hsv(8), ...
        'FaceColor',RB(j).color )
    
    patch( 'Vertices',RBc(j).conf, ...
        'Faces',faces_matrix,...
        'FaceVertexCData',hsv(8), ...
        'FaceColor',RBc(j).color )
    
end

c = [0.749019622802734 0 0.749019622802734] ;
id1 = [1 2 3]; 
id2 = [1 2 4]; 
id3 = [2 3 4]; 
id4 = [3 1 4];

for j = 1 : nbp-1
    
    
    patch(Pw(id1, 1,j), Pw(id1, 2,j), Pw(id1, 3,j), c ) 
    patch(Pw(id2, 1,j), Pw(id2, 2,j), Pw(id2, 3,j), c )
    patch(Pw(id3, 1,j), Pw(id3, 2,j), Pw(id3, 3,j),'m')
    patch(Pw(id4, 1,j), Pw(id4, 2,j), Pw(id4, 3,j), c )
    
    patch(Pc(id1, 1,j), Pc(id1, 2,j), Pc(id1, 3,j), c )
    patch(Pc(id2, 1,j), Pc(id2, 2,j), Pc(id2, 3,j), c )
    patch(Pc(id3, 1,j), Pc(id3, 2,j), Pc(id3, 3,j),'m')
    patch(Pc(id4, 1,j), Pc(id4, 2,j), Pc(id4, 3,j), c )
end



end

