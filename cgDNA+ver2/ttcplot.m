function Data = ttcplot(Data)

nbr_Seq = length(Data);

j = 1 ; 


for i = 1:nbr_Seq
    
    c = rand(1,3) ;
    
    if isfield(Data,'ttc')~=1 && isfield(Data,'ttc_intr')~=1
        
        
        Data = cgDNApttc(Data,10000,0);
        
    end
    
    ttc = Data(i).ttc' ;
    ttc_intr = Data(i).ttc_intr' ;
    nbp = length(ttc);
    
    plot(1:nbp, log(ttc), '-', 'color', c, 'LineWidth', 2 );
    
    hold on
    plot(1:nbp, log(ttc) - log(ttc_intr), '--','color', c, 'LineWidth', 2 );
    
    lg{j} = Data(i).seqlabel ; 
    lg{j+1} = [ Data(i).seqlabel ' fact'  ] ; 
    
    j = j + 2; 
end

xlabel('basepair level')

legend(lg,'FontSize',12)

hold off

end