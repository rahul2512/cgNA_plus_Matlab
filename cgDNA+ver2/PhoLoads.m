function Data = PhoLoads(input,sequence)

Path_to_Libraries1 = '/Users/patelli/Documents/cgDNAMechanics' ;
addpath(genpath(Path_to_Libraries1),genpath('./')) ;

if isstring(input)
    
    shape = LoadFrames( input , sequence );
    
else
    shape = input;
end

Data = cgDNAp(sequence);
dE = Data.stiff * (shape - Data.groundstate) ;

Data(2) = Data(1) ;
Data(2).seqlabel = 'Loads';
Data(2).groundstate = shape;

[~,~,dEetapW,dEwpW,~,~,dEetapC,dEwpC] = vector2shapes(dE);

[~,~,etapW,wpW,~,~,etapC,wpC] = vector2shapes(shape);

n = length(sequence);

uscale = 5;

%% Compute lam_W
for m = 1:n-1
    
    Bpm_R = Cayleyu2Q(etapW(m,:)',uscale);
    Bpm_r = wpW(m,:)';
    
    Ad_Bpm = Admat(Bpm_R,Bpm_r);
    
    L = [ Pmat(etapW(m,:),uscale) , zeros(3,3) ; zeros(3,3) eye(3) ] ;
    
    lam_W(:,m) = - Ad_Bpm'\(L'\[dEetapW(m,:)' ; dEwpW(m,:)']) ;
    
end

%% Compute lam_C
for m = 1:n-1
    
    Bpm_R = Cayleyu2Q(etapC(m,:)',uscale);
    Bpm_r = wpC(m,:)';
    
    Ad_Bpm = Admat(Bpm_R,Bpm_r);
    
    L = [ Pmat(etapC(m,:),uscale) , zeros(3,3) ; zeros(3,3) eye(3) ] ;
    
    lam_C(:,m) = - Ad_Bpm'\(L'\[dEetapC(m,:)' ; dEwpC(m,:)']) ;
    
end

Data(2).lam_W = lam_W;
Data(2).lam_C = lam_C;

rmpath(genpath(Path_to_Libraries1))

end