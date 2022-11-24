close all
n = 5;

load 'olig_Palin_141_3us_symmetrized_ver2.mat'

seq = olig(n).seq;

Data = cgDNAp(seq);

if 0
    
    Data(2) = Data(1);
    Data(2).seqlabel = 'MD';
    Data(2).groundstate = olig(n).shape_sym;
    Data(2).stiff = olig(n).stiff_me_sym;
    
    cgDNAp2dplot(Data,'RotDeg',0)
    close(figure(1))
    
elseif 1
    
    Data.groundstate = abs(Data.groundstate-olig(n).shape_sym);
    Data.stiff = olig(n).stiff_me_sym;
    
    Data(2) = Data(1);
    
    addpath(genpath('/Users/patelli/Documents/MATLAB/cgDNA'));
    ps = load('cgDNAparamset4.maxent.mat');

    [shapes, ~] = constructGaussian(seq, ps);
    
    shapes = reshape([ shapes ; nan(6,1) ], [12 24 ] )';
    
    shapes = [shapes nan(24,12)];
    shapes = shapes(: , [ 1:6, 13:18 ,7:12 , 19:24 ] )';
    shapes = shapes(:);
    shapes = shapes(1:end-18) ;
    Data(2).groundstate = abs(shapes-olig(n).shape_sym);
    Data(2).seqlabel = 'cgDNA';
    Data(1).seqlabel = 'cgDNA+';
   
    rmpath('/Users/patelli/Documents/MATLAB/cgDNA');
    
    cgDNAp2dplot(Data,'RotDeg',0)
    close(figure(2))
    
end