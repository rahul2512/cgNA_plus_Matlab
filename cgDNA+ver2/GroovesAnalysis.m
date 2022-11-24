function GroovesAnalysis(Data)

figure 
addpath('./Utilities/')

nbr_Seq = length(Data) ;

if nbr_Seq == 1
    
    D = getDistanceMatrix(Data.groundstate);
    plotDistanceMatric(D,Data) ;
        
elseif nbr_Seq > 1 && nbr_Seq <=4
    
    for i = 1:nbr_Seq
        D = getDistanceMatrix(Data(i).groundstate);
        subplot(ceil(nbr_Seq/2),2,i)
        plotDistanceMatric(D,Data(i),'subplot') ;
    end
    
    
else
    error('This function works for a maxium of four sequences.');
    
    
end

end

function D = getDistanceMatrix(shapes)

bs_level = frames(shapes);
pW = [bs_level.rpw];
pC = [bs_level.rpc];


f1=fnplt(cscvn(pW));
f2=fnplt(cscvn(pC));

D = pdist2(f1',f2');

end
function plotDistanceMatric(D,myData,~)

if nargin > 3 
    ft_size = 12 ; 
else 
    ft_size = 15 ; 
end

[C,h]=contour(D,'LevelStep',1,'TextStep',2,'ShowText','on','Fill','on','LineColor','k');
clabel(C,h,'FontSize',30,'FontWeight','bold','Color',[0 0 0]);
c = [ linspace(0,1,64)' linspace(0,1,64)' ones(64,1)  ] ;
colormap(c);
title(myData.sequence,'FontSize',30)
xlabel('Reading strand','FontSize',20)
ylabel('Complementary strand','FontSize',20)
axis square

end
function c = mycolormap()

c = [         0         0    0.5625
         0         0    0.5938
         0         0    0.6250
         0         0    0.6562
         0         0    0.6875
         0         0    0.7188
         0         0    0.7500
         0         0    0.7812
         0         0    0.8125
         0         0    0.8438
         0         0    0.8750
         0         0    0.9062
         0         0    0.9375
         0         0    0.9688
         0         0    1.0000
    0.0471    0.0471    0.9882
    0.0941    0.0941    0.9765
    0.1412    0.1412    0.9647
    0.1882    0.1882    0.9529
    0.2353    0.2353    0.9412
    0.2824    0.2824    0.9294
    0.3294    0.3294    0.9176
    0.3765    0.3765    0.9059
    0.4235    0.4235    0.8941
    0.4706    0.4706    0.8824
    0.5176    0.5176    0.8706
    0.5647    0.5647    0.8588
    0.6118    0.6118    0.8471
    0.6588    0.6588    0.8353
    0.7059    0.7059    0.8235
    0.7529    0.7529    0.8118
    0.8000    0.8000    0.8000
    0.8133    0.7467    0.7467
    0.8267    0.6933    0.6933
    0.8400    0.6400    0.6400
    0.8533    0.5867    0.5867
    0.8667    0.5333    0.5333
    0.8800    0.4800    0.4800
    0.8933    0.4267    0.4267
    0.9067    0.3733    0.3733
    0.9200    0.3200    0.3200
    0.9333    0.2667    0.2667
    0.9467    0.2133    0.2133
    0.9600    0.1600    0.1600
    0.9733    0.1067    0.1067
    0.9867    0.0533    0.0533
    1.0000         0         0
    0.9706         0         0
    0.9412         0         0
    0.9118         0         0
    0.8824         0         0
    0.8529         0         0
    0.8235         0         0
    0.7941         0         0
    0.7647         0         0
    0.7353         0         0
    0.7059         0         0
    0.6765         0         0
    0.6471         0         0
    0.6176         0         0
    0.5882         0         0
    0.5588         0         0
    0.5294         0         0
    0.5000         0         0] ;


end

