Path_to_Libraries1 = '../../../CorePart' ;
Path_to_Libraries2 = '../../Library' ;
addpath(genpath(Path_to_Libraries1),genpath(Path_to_Libraries2),genpath('./')) ;
load('./Data/PDBatoms.mat')

[ConstitutiveParam, basepair ] =prep_ConstitutiveParam(basepair_fit) ;

clear coord r_factor basepair_fit

n = 10 ;

c = c(n,:,2)' ; f = f(n,:,2)' ;

shape = ChainCoord(basepair);

[eta, w, u, v] = DecompCoord(shape) ;

dE=gradEnergyBichain(ConstitutiveParam,shape);
[dEeta, dEw, dEu, dEv] = DecompCoord(dE) ;

%% Compute Zeta_2
Pn = junctConf(eta(n,:),w(n,:));

Pn_rot = Pn(1:3,1:3);
Pn_half_rot = sqrtm(Pn(:,1:3));

mn = ( Pn_half_rot + eye(3) )\eye(3) ;
  
Ad_P_h = Admat(Pn_half_rot',-0.5*w(n,:));

X = 0.5*(vect2mat(w(n,:))*( Pn_half_rot' + eye(3)) )*mn ;

A = [-Pn_half_rot'*mn , zeros(3,3) ; X , -0.5*Pn_half_rot'] ;
Ad_P = Admat(Pn);
 
Fn = Ad_P_h + A ; 

Zeta_2 = 0 ; 
if n<10
    dExn = [dEu(n,:) , dEv(n,:)]' ;
    L_xn = Linvmat([u(n,:) , v(n,:)],5);
    Zeta_2 = L_xn'*dExn ;
end

%% Compute Zeta_1
Zeta_1 = 0 ;
if n > 1
    anm1 = junctConf(u(n-1,:),v(n-1,:));
    Ad_anm1 = Admat(anm1);
    dExnm1 = [dEu(n-1,:) , dEv(n-1,:)]' ;
    L_xnm1  = Linvmat([u(n-1,:) , v(n-1,:)],5);
    Zeta_1 = Ad_anm1'*L_xnm1'*dExnm1 ;
end

%% Compute mu_P
dEyn = [dEeta(n,:) , dEw(n,:)]' ; 
L_yn_Ene = Linvmat([eta(n,:) , w(n,:)],5)' * dEyn;

Ad_Pn = Admat(Pn);
mu_P = L_yn_Ene ;

stress_plus = Fn'*(-Zeta_2 + Zeta_1) - mu_P;
stress_plus = -stress_plus';
Load = [ c ; f ]' ;

clearvars -except stress_plus Load
