function shape = LoadFrames(filename, sequence)

nbp = length(sequence);

[ base , pho ] = importFrames(filename);

baseW = base(1:nbp,:);
baseC = base(nbp+1:2*nbp,:);

phoW = pho(1:(nbp-1),:);
phoC = pho(nbp:2*(nbp-1),:);

for i = 1:nbp
    basepair(i).Rw = reshape(baseW(i,1:9), [3 3]);
    basepair(i).rw = baseW(i,10:12)';
    
    Rc = reshape(baseC(i,1:9), [3 3])*diag([1 -1 -1]);
    
    basepair(i).Rc = Rc;
    basepair(i).rc = baseC(i,10:12)';
    
    basepair(i).R =  Rc*sqrtm(Rc'*basepair(i).Rw);
    basepair(i).r = 0.5*(basepair(i).rw + basepair(i).rc);
    
    basepair(i).seq = sequence(i);
    
end

for i = 1:nbp-1
    basepair(i+1).Rpw = reshape(phoW(i,1:9), [3 3]);
    basepair(i+1).rpw = phoW(i,10:12)';
    basepair(i).Rpc = reshape(phoC(i,1:9), [3 3]);
    basepair(i).rpc = phoC(i,10:12)';
end

eta = zeros(nbp,3);
w = zeros(nbp,3);
u = zeros(nbp-1,3);
v = zeros(nbp-1,3);

etapW = zeros(nbp-1,3);
wpW = zeros(nbp-1,3);

etapC = zeros(nbp-1,3);
wpC = zeros(nbp-1,3);

uscale = 5 ;

for i = 1:(nbp - 1)
    
    [Q,q]=SEmult(SEinv([basepair(i).R basepair(i).r]),[basepair(i+1).R basepair(i+1).r]);
    [u(i,:),v(i,:)]=junctCoord(Q,q,uscale);
    
end

for j = 1:nbp
    
    [Q,q]=SEmult(SEinv([basepair(j).Rc basepair(j).rc]),[basepair(j).Rw basepair(j).rw]);
    [eta(j,:),w(j,:)] = junctCoord(Q,q,uscale) ;
    
end

for k = 1:nbp-1
    j = k+1;
    [Qw,qw]=SEmult(SEinv([basepair(j).Rw basepair(j).rw]),[basepair(j).Rpw basepair(j).rpw]);
    etapW(k,:) = uscale*(2/(trace(Qw)+1))*mat2vect(Qw);
    wpW(k,:)   = qw ;
    
    [Qc,qc]=SEmult(SEinv([basepair(k).Rc*diag([1 -1 -1]) basepair(k).rc]),[basepair(k).Rpc basepair(k).rpc]);
    etapC(k,:) = uscale*(2/(trace(Qc)+1))*mat2vect(Qc);
    wpC(k,:)   = qc ;
    
end

N = 24*nbp - 18;

rotat = [ eta ; etapC ; u ; etapW ];
trasl = [   w ; wpC ; v ; wpW ];
q = [ rotat, trasl ];
indices = repmat( 1:nbp-1 , [4 1] ) + repmat( [ 0 ; nbp ; 2*(nbp-1)+1 ; 3*(nbp-1)+1] , [1 nbp-1] ) ;
indices = [indices(:); nbp];

shape = reshape( q(indices,:)', 1, N )';

end

function a=mat2vect(A)
%
% Compute the triple a associated with
% the skew part of the matrix A

a=zeros(1,3);

a(1,1)=A(3,2)-A(2,3);
a(1,2)=A(1,3)-A(3,1);
a(1,3)=A(2,1)-A(1,2);

end

