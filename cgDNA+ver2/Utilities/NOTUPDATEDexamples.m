%---------------------------------------------------------
% Welcome to the cgDNA package!
% 
% Run this example script to discover all the possible cgDNA functionalities
% (coordinate plots, 3D viewer, persistence length computations, etc).
% 
% The main function of this program is cgDNA.m. It can be run with only one sequence
% (EXAMPLE 1 below) or multiple sequences (EXAMPLE 2 below) as an input.
% 
% Have fun testing this script with different sequences as input, and/or use
% it as a template for your own code/computations.
%---------------------------------------------------------


%---------------------------------------------------------
%% EXAMPLE 1: single sequence
%---------------------------------------------------------
%----------------------------------------------------------
% ... Pick a sequence ...

sequence = 'CGCGATATCGCG';

fprintf('Working on sequence %s \n',sequence);

% ... or read the sequence from a plain text with the sequence 
% or from a FASTA format file

%sequence = parseFasta('sequence.txt');

% ... or build a repeat using the poly function (can be nested)...
sequence_repeat = poly(sequence, 10);

%----------------------------------------------------------
% ... Construct cgDNA groundstate vector and stiffness matrix ...

fprintf('Building cgDNA data (groundstate vector and stiffness matrix)... \n');

myData =  cgDNAp(sequence);

% ... You can specify a cgDNA parameter set (default is 4) and label your sequence:
% ... (see cgDNA.m for details)

myData2 = cgDNAp(sequence_repeat, 1, 'myRepeat');

%----------------------------------------------------------
% ... plot the groundstate coordinates ...

fprintf('plotting groundstate coordinates... \n')
 cgDNAp2dplot(myData);
%----------------------------------------------------------
% ... 3D view of the groundstate ... 

fprintf('plotting groundstate 3D view... \n')
 cgDNAp3dplot(myData);
%----------------------------------------------------------
% ... Monte Carlo sampling of the cgDNA probability density ...

 N_samp = 100; %pick of number of samples 
 fprintf(['Sampling ' num2str(N_samp) ' oligomer conformations of the sequence %s: \n'])
 fprintf( [myData2.sequence ' \n' ])
 
 myData2 = cgDNAp_MonteCarlo(myData2,N_samp);
 
% ... from which persistence lengths are estimated ... 
 
 fprintf(['apparent persistence length: ' num2str(myData2.apparent_pl) 'bp \n' ])
 fprintf(['dynamic persistence length: ' num2str(myData2.dynamic_pl) 'bp \n' ])

%----------------------------------------------------------  
% ... Construct PDB file of atomic coordinates ...
% for each base and save results in a file.
% NOT UP-TO-DATE 
% PDBOutputFile = 'base_atoms.pdb';
% fprintf('Saving coordinates to PDB file <%s>... \n', PDBOutputFile);
% makePDB(myData, PDBOutputFile);

%----------------------------------------------------------
%---------------------------------------------------------- 



% %---------------------------------------------------------
% %% EXAMPLE 2: multiple sequences
% %---------------------------------------------------------
% %---------------------------------------------------------
% % ... Multiple sequences can be input as a cell array ...
% % 
% %  sequences = {  'CGCGAATTCGCG';
% %                 'CGCGAGTT';
% %                  poly('A',10)
% %                 };
%             
% % ... with parameter sets and labels ...       
%  sequences = {  'CGCGAATTCGCG', 4, 'a_sequence';
%                 'CGCGAATT',3, 'a_shorter_one';
%                  poly('A',10), 4, 'PolyA'
%                 };    
%  
% 
% % ... or , if you have access to the MATLAB Bioinformatics toolbox,
% % read the sequences from a FASTA format file 
% % (parameter sets must be added manually via a cell array)
% 
% % fasta_struct = fastaread('sequence_mult.fa');
% % sequences = [{fasta_struct.Sequence}', {4;3;4}, {fasta_struct.Header}'];
%  
%          
% %----------------------------------------------------------
% % ... Construct cgDNA groundstate vector and stiffness matrix ..
% 
%  fprintf('Building cgDNA data (groundstate vector and stiffness matrix)... \n');
%  myData_mult =  cgDNA(sequences);
%  
% 
% %----------------------------------------------------------
% % ... plot all the groundstate coordinates ...
% 
%  fprintf('plotting groundstate coordinates... \n')
%  cgDNA2dplot(myData_mult);
% %------------------s----------------------------------------
% % ... 3D view of all the groundstates ... 
% 
%  fprintf('plotting groundstate 3D view... \n')
%  cgDNA3dplot(myData_mult);
% %----------------------------------------------------------
% % ... Monte Carlo sampling of the cgDNA probability densities ...
% % 
% 
%  fprintf(['Sampling ' num2str(N_samp) ' oligomer configurations of the sequences \n'])
%  N_samp = 100; %pick of number of samples 
%  myData_mult = cgDNA_MonteCarlo(myData_mult,N_samp);
%  
% % ... from which persistence lengths are estimated ... 
%   for i= 1:length(myData_mult)
%       fprintf(['sequence: ' myData_mult(i).sequence ' \n' ])
%       fprintf(['apparent persistence length: ' num2str(myData_mult(i).apparent_pl) 'bp \n' ])
%       fprintf(['dynamic persistence length: ' num2str(myData_mult(i).dynamic_pl) 'bp \n' ])
%   end
%  
%   
% %----------------------------------------------------------  
% % ... Construct PDB file of atomic coordinates ...
% % for each base and save results in a file.
% PDBOutputFile = 'base_atoms.pdb';
% fprintf('Saving coordinates to PDB file <%s>... \n', PDBOutputFile);
% makePDB(myData_mult, PDBOutputFile);
% 
% %---------------------------------------------------------

