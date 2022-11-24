function Data = SugarPseudorot(Data)

SugarAtomsName={'O4''','C1''','C2''','C3''','C4'''};
List = ListOfDihedralAngle();

W = [Data.AtomsCoord.W];
C = [Data.AtomsCoord.C];

nbp = Data.nbp ;

for i = 1:nbp
  
    [~, idx]=ismember(W(i).names, SugarAtomsName);
    
    AtomsCoord = W(i).coord(:,idx~=0);
    [~,trf] = sort(idx(idx~=0), 'ascend');
    AtomsCoord = AtomsCoord(:,trf);
    
    A = 0 ; 
    B = 0 ;
    for ang = 1:5
        phi = getDihedralAngle(AtomsCoord(:,List(ang,:)));
        
        A = A + phi * cos(4/5*pi*(ang-1));
        B = B + phi * sin(4/5*pi*(ang-1));
        
    end
    A =  2/5*A;
    B = -2/5*B;
    
    % TODO : Change according to Wolfram Saenger Page 20
    
    phase = atan2(B,A)*180/pi;
    
    if phase < 0 
        phase = phase + 360;
    end
    
    if phase < 0 
        phase = phase + 360;
    end
    
    if phase < 36
		pucker = 'C3''-endo';
    elseif phase < 72
        pucker = 'C4''-exo';
	elseif phase < 108
		pucker = 'O4''-endo';
	elseif phase < 144
		pucker = 'C1''-exo';
	elseif phase < 180
		pucker = 'C2''-endo';
	elseif phase < 216
		pucker = 'C3''-exo';
	elseif phase < 252
		pucker = 'C4''-endo';
	elseif phase < 288
		pucker = 'O4''-exo';
	elseif phase < 324
		pucker = 'C1''-endo';
	elseif phase < 360
		pucker = 'C2''-exo';
    end
    
    
    Data.AtomsCoord(i).W.phase = phase ;
    Data.AtomsCoord(i).W.amp = sqrt(A^2+B^2); 
    Data.AtomsCoord(i).W.pucker = pucker ;
    
end

for i = 1:nbp
  
    [~, idx]=ismember(C(i).names, SugarAtomsName);
    
    AtomsCoord = C(i).coord(:,idx~=0);
    [~,trf] = sort(idx(idx~=0), 'ascend');
    AtomsCoord = AtomsCoord(:,trf);
    
    A = 0 ; 
    B = 0 ;
    for ang = 1:5
        phi = getDihedralAngle(AtomsCoord(:,List(ang,:)));
  
        A = A + phi * cos(4/5*pi*(ang-1));
        B = B + phi * sin(4/5*pi*(ang-1));
        
    end
    A =  2/5*A;
    B = -2/5*B;
    
    phase = atan2(B,A)*180/pi; 
    
    if phase < 0 
        phase = phase + 360;
    end
    
    if phase < 36
		pucker = 'C3''-endo';
    elseif phase < 72
        pucker = 'C4''-exo';
	elseif phase < 108
		pucker = 'O4''-endo';
	elseif phase < 144
		pucker = 'C1''-exo';
	elseif phase < 180
		pucker = 'C2''-endo';
	elseif phase < 216
		pucker = 'C3''-exo';
	elseif phase < 252
		pucker = 'C4''-endo';
	elseif phase < 288
		pucker = 'O4''-exo';
	elseif phase < 324
		pucker = 'C1''-endo';
	elseif phase < 360
		pucker = 'C2''-exo';
    end
    
    Data.AtomsCoord(i).C.phase = phase ;
    Data.AtomsCoord(i).C.amp = sqrt(A^2+B^2); 
    Data.AtomsCoord(i).C.pucker = pucker ;
    
    
end


end

function List = ListOfDihedralAngle()

List = nan(5,4); 
List(1,:) = [2 3 4 5];
List(2,:) = [3 4 5 1];
List(3,:) = [4 5 1 2];
List(4,:) = [5 1 2 3];
List(5,:) = [1 2 3 4];

end

function phi = getDihedralAngle(R)

r1 = R(:,1) ; r2 =R(:,2) ; r3 = R(:,3); r4 = R(:,4) ;

b1 = r2 - r1 ; b2 = r3-r2 ; b3 = r4 - r3 ;

n1 = cross(b1,b2); n2 = cross(b2,b3);

m1 = cross(n1,b2/norm(b2));

x = n1'*n2 ; y = m1'*n2 ;

phi = -180/pi*atan2( y,x ) ;


end

