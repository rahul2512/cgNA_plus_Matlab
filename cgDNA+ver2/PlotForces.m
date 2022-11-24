function PlotForces(Data)

lam_W = Data(2).lam_W;
lam_C = Data(2).lam_C; 

subplot(3,2,1)
plot(lam_W(1:3,:)') ; 
legend({'1','2','3'})
title('Couple Watson')

subplot(3,2,2)
plot(lam_W(4:6,:)') ; 
legend({'1','2','3'})
title('Force Watson')

subplot(3,2,3)
plot(lam_C(1:3,:)') ;
legend({'1','2','3'})
title('Couple Crick')

subplot(3,2,4)
plot(lam_C(4:6,:)') ;
legend({'1','2','3'})
title('Force Crick')

subplot(3,2,5)
plot( diag(sqrt(lam_W(1:3,:)'*lam_W(1:3,:))), '-k' ) ; hold on 
plot( diag(sqrt(lam_C(1:3,:)'*lam_C(1:3,:))), '--k' ) ; hold off
legend({'Watson','Crick'})
title('Norm Couple')

subplot(3,2,6)
plot( diag(sqrt(lam_W(4:6,:)'*lam_W(4:6,:))), '-k') ; hold on
plot( diag(sqrt(lam_C(4:6,:)'*lam_C(4:6,:))), '--k') ; hold off
legend({'Watson','Crick'})
title('Norm Force')


end

