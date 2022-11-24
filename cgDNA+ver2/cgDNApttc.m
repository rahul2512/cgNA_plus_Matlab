function Data = cgDNApttc(Data,NbrSample,Base2Drop,varargin)

addpath('./Utilities')
addpath(genpath('./cgDNApmc'))

% Check varargin for optional input
p = CheckInput(varargin);

NS = num2str(NbrSample);
B2D =num2str(Base2Drop);

nbr_Seq = length(Data);

home = pwd;

PathToFile  = './cgDNApmc_files/cgDNApmc';

for i = 1:nbr_Seq
    
    RunName = [ Data(i).seqlabel '_' num2str(randi([0 999999] )) ] ;
    
    cd(PathToFile)
    
    cgDNApmc('t0',RunName,Data(i).sequence,'cgDNA+ps1.txt',p.GenOrAcc,NS,B2D,p.J) ;
    
    [bp, ttc, ttc_intr] = importOutput(RunName) ;
    
    [pl, pl_fact, ~] = Compute_pl(length(bp), ttc, ttc_intr) ;
    
    Data(i).ttc = ttc ; 
    Data(i).ttc_intr = ttc_intr ;
    Data(i).apparent_pl = pl ;
    Data(i).dynamic_pl  = pl_fact;
    
    fprintf('Computation for Seq %d done \n', i);
    
    cd(home)
    RunName_Old{i} = RunName ; 
    
end

RemoveOutput(RunName_Old,PathToFile);

end

function p = CheckInput(inputarg)

p = inputParser ;

validGen = @(x) ischar(x) && contains('ga',x);
validJ = @(x) ischar(x) && contains('ny',x);

addOptional(p,'GenOrAcc','g',validGen) ;
addOptional(p,'J','n',validJ) ;

parse(p,inputarg{:});

p=p.Results;


end
function [pl, pl_fact, tmp_pl] = Compute_pl(nbp,dot_d3,dot_d3_intr)

warning('off')

x = [ zeros(nbp,1) , (1:nbp)' ];
y = log(dot_d3) ;

tmp_pl = x \ y ;
pl = -1/tmp_pl(2) ; 

y = log(dot_d3) - log(dot_d3_intr) ;

tmp_pl = x \ y ;
pl_fact = -1/tmp_pl(2) ; 
end
function RemoveOutput(RunName_Old,PathToFile)

for i =1:length(RunName_Old)
   
    cmd = [ 'rm ' PathToFile '/' RunName_Old{i} '*' ] ; 
    system(cmd) ;
    
end

end

