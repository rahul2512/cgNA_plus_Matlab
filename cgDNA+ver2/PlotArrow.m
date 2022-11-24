function PlotArrow(Data)

addpath(genpath('./Utilities'))

cgDNAp3dplot(Data)
hold on
bp_level = frames(Data.groundstate);

nbp = length(bp_level) ;
arr_length=5;
map = colormap('copper') ;

for k = 1:2
    
    range = 1+3*(k-1):3+3*(k-1);
    
    D1 = Data.lam_W(range,:) ;
    D2 = Data.lam_C(range,:) ;
    
    M1 = vecnorm(D1) ;
    M2 = vecnorm(D2) ;
    
    M = max( max(M1), max(M2)  )  ;
    
    ax1(k) = subplot(1, 2, k);
    cgDNApviewer(Data.groundstate, Data.sequence,1) ;
    for i = 1:nbp-1
        
        h = D1(:,i)/M1(i)*arr_length ;
        O = bp_level(i+1).rpw ;
        P = O + h ;
        mArrow3(O,P,'tipWidth',0.3,'color',map(ceil(M1(i)/M*length(map)),:) ) ;
        
    end
    
    for i = 1:nbp-1
        
        h = D2(:,i)/M2(i)*arr_length ;
        O = bp_level(i).rpc;
        P = O + h ;
        mArrow3(O,P,'tipWidth',0.3,'color',map(ceil(M2(i)/M*length(map)),:) ) ;
        
    end
    
    axis equal
    
    colorbar
    caxis([0 M])

end

Link = linkprop([ax1(1), ax1(2)],{'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);

hold off

end