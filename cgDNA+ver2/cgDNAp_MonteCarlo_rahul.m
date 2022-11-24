function Data = cgDNAp_MonteCarlo(Data,nbr_samples)
% ======================================================================== 
% cgDNAp_MonteCarlo(Data,nbr_samples)
% 
% This function draws a given number of configurations from the cgDNA 
% probability  density function for the sequence(s) in the given Data
% structure.
%
% Apparent and dynamic (see ref below for definitions) persistence lengths are computed from the 
% sampling for each sequence, and are added as fields in the structure
% Data.
%
% This is a MATLAB version of the much more efficient cgDNAmc code 
% http://lcvmwww.epfl.ch/software/cgDNAmc/doc/index.html
% Jonathan S. Mitchell, Jaroslaw Glowacki, Alexandre E. Grandchamp, Robert S. Manning, and John H. Maddocks, 
% Sequence-dependent persistence lengths of DNA 
% J. Chem. Theory Comput.
%
% 
%
%
% If only one sequence is provided as input, this function also plots all the  
% base-pairs positions of all sampled configuration in blue and the 
% ground-state in red.
% 
% Input : - Data           Data structure as output by the cgDNA.m function
%         - nbr_samples    Number of wanted sampled configuration. We
%                          advise to not use above 500 samples. 
% Output : - Data          Input Data structure, to which fields and values
%                          for apparent and dynamic persistence length were added
%                          for each sequence
%

addpath('./Utilities');

L_data = length(Data);

for j =1:L_data
    
    % Load the cgDNA Data
    params = load(['./ParameterSets/' Data(j).paramset '.mat']) ;
    sequence = Data(j).sequence;

    nbp = Data(j).nbp;

    w = Data(j).groundstate;
    K = Data(j).stiff;
    dim = length(w) ;

    % Compute the Cholesky factorisation of K
    L = chol(K) ;

    % Sample from the multivariate standard normal distribution
    MC = mvnrnd(zeros(dim,1),eye(dim),nbr_samples) ;

    if L_data == 1
        figure
    end
    
    ttc = zeros(nbp,1);
    % Loop over the samples
    rand = zeros(dim,nbr_samples);
    for i = 1:nbr_samples

      % Solve the linear system to get cgDNA coordinates
      x = L \ MC(i,:)' + w ; 
      rand(:,i) = x ;  
      % Reconstruct the frames for the coordinate x
      bp_lv = frames(x) ;

      XYZ = zeros(nbp,3);

      % Loop over the basepair
      for k = 1:nbp

        % Store all the xyzc oordinates of the base-pair positions
        XYZ(k,:) = bp_lv(k).r ;
        ttc(k,1) = ttc(k,1) + bp_lv(k).R(3,3)/nbr_samples;

      end

      if L_data == 1
        plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'-b'), hold on ;
      end

    end

     % Reconstruct the frames for the ground state 
    bp_lv = frames(w) ;

    % Loop over the basepair

    ttc_int = zeros(nbp,1);
    for k = 1:length(sequence)

     % Store all the xyzc oordinates of the base-pair positions
      XYZ(k,:) = bp_lv(k).r ;
      ttc_int(k,1)= bp_lv(k).R(3,3);

    end

    %Compute persistence length
    [l_p{j}, l_d{j}, ~] = Compute_pl(nbp,ttc,ttc_int);


    if L_data == 1
        plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'-xr')

        axis equal 
        grid on 

        hold off 
    end
    
end
    
% Add persistence lengths value to Data structure
    [Data.apparent_pl] = l_p{:};
    [Data.dynamic_pl] = l_d{:};
    [Data.config] = rand;
   

end


% Fit persistence length to tangent-tangent correlation
function [pl, pl_fact, tmp_pl] = Compute_pl(nbp,dot_d3,dot_d3_intr)

warning('off')

x = [ zeros(nbp,1) , [1:nbp]' ];
y = log(dot_d3) ;

tmp_pl = x \ y ;
pl = -1/tmp_pl(2) ; 

y = log(dot_d3) - log(dot_d3_intr) ;

tmp_pl = x \ y ;
pl_fact = -1/tmp_pl(2) ; 


end