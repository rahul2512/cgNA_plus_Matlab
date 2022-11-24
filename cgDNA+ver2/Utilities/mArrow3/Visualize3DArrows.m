function Visualize3DArrows( basepair ,ConstitutiveParam )



[c,f]= LoadsDoubleChain(basepair, ConstitutiveParam) ;


plotArrow(f,basepair)
plotArrow(c,basepair)


end

function plotArrow(data,basepair)

figure

nbp = length(basepair) ;
arr_length=5;
map = colormap(cool) ;

D1 = squeeze(data(:,1,:)) ;
D2 = squeeze(data(:,2,:)) ;

M1 = vecnorm(D1) ;
M2 = vecnorm(D2) ; 

M = max( max(M1), max(M2)  )  ;

for i = 1:nbp
  
  h = D1(:,i)/M1(i)*arr_length ;
  O = (basepair(i).rw + basepair(i).Rw*[ 0 2.5 0 ]' ) ;
  P = O + h ;
  mArrow3(O,P,'stemWidth',0.15,'tipWidth',0.3,'color',map(ceil(M1(i)/M*length(map)),:) ) ;
  
  hold on
  
  h = D2(:,i)/M2(i)*arr_length ;
  O = (basepair(i).rc - basepair(i).Rc*[ 0 2.5 0 ]'  ) ;
  P = O + h ;
  mArrow3(O,P,'stemWidth',0.15,'tipWidth',0.3,'color',map(ceil(M2(i)/M*length(map)),:) ) ;
end

axis equal

colorbar
caxis([0 M])

ChainViewer(basepair)
hold off


end



