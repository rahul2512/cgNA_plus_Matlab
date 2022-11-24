function Data = BI_BIIstates(Data,pos)
    
Epsilon_1 = {'C4''','C3''','O3'''}; % P_(i+1)
Epsilon_2 = 'P  ';
Zeta_1    = {'C3''','O3'''};
Zeta_2    = {'P  ', 'O5'''}; % P_(i+1) O5'_(i+1)

W = [Data.AtomsCoord.W];
C = [Data.AtomsCoord.C];

nbp = Data.nbp ;

if nargin == 1 
    pos = 1:nbp-1;
end

for i = pos
  
    
    % Get Epsilon torsional angle
    [~, idx]=ismember(W(i).names, Epsilon_1);
    AtomsCoord = W(i).coord(:,idx~=0);
    [~,trf] = sort(idx(idx~=0), 'ascend');
    AtomsCoord = AtomsCoord(:,trf);
    [~, idx]=ismember(W(i+1).names, Epsilon_2);
    AtomsCoord = [ AtomsCoord W(i+1).coord(:,idx~=0) ];
    
    phi_eps = getDihedralAngle(AtomsCoord);

    % Get Zeta torsional angle
    [~, idx]=ismember(W(i).names, Zeta_1);
    AtomsCoord = W(i).coord(:,idx~=0);
    [~,trf] = sort(idx(idx~=0), 'ascend');
    AtomsCoord = AtomsCoord(:,trf);
    [~, idx]=ismember(W(i+1).names, Zeta_2);
    AtomsCoord = [ AtomsCoord W(i+1).coord(:,idx~=0) ];
    
    phi_zeta = getDihedralAngle(AtomsCoord);
    
    if phi_eps - phi_zeta < 0
        
        state = 'BI';
        
    elseif phi_eps - phi_zeta > 0
        
        state = 'BII';
        
    end
    
    Data.AtomsCoord(i).W.epsilon = phi_eps;
    Data.AtomsCoord(i).W.zeta    = phi_zeta;
	Data.AtomsCoord(i).W.state   = state;   
        
end

for i = pos
  
    % Get Epsilon torsional angle
    [~, idx]=ismember(C(i).names, Epsilon_1);
    AtomsCoord = C(i).coord(:,idx~=0);
    [~,trf] = sort(idx(idx~=0), 'ascend');
    AtomsCoord = AtomsCoord(:,trf);
    [~, idx]=ismember(C(i+1).names, Epsilon_2);
    AtomsCoord = [ AtomsCoord C(i+1).coord(:,idx~=0) ];
    
    phi_eps = getDihedralAngle(AtomsCoord);

    % Get Zeta torsional angle
    [~, idx]=ismember(C(i).names, Zeta_1);
    AtomsCoord = C(i).coord(:,idx~=0);
    [~,trf] = sort(idx(idx~=0), 'ascend');
    AtomsCoord = AtomsCoord(:,trf);
    [~, idx]=ismember(C(i+1).names, Zeta_2);
    AtomsCoord = [ AtomsCoord C(i+1).coord(:,idx~=0) ];
    
    phi_zeta = getDihedralAngle(AtomsCoord);
    
    if phi_eps - phi_zeta < 0
        
        state = 'BI';
        
    elseif phi_eps - phi_zeta > 0
        
        state = 'BII';
        
    end
    
    Data.AtomsCoord(i).C.epsilon = phi_eps;
    Data.AtomsCoord(i).C.zeta    = phi_zeta;
	Data.AtomsCoord(i).C.state   = state;
    
end

end

function phi = getDihedralAngle(R)

r1 = R(:,1) ; r2 =R(:,2) ; r3 = R(:,3); r4 = R(:,4) ;

b1 = r2 - r1 ; b2 = r3-r2 ; b3 = r4 - r3 ;

n1 = cross(b1,b2); n2 = cross(b2,b3);

m1 = cross(n1,b2/norm(b2));

x = n1'*n2 ; y = m1'*n2 ;

phi = -180/pi*atan2( y,x ) ;


end

