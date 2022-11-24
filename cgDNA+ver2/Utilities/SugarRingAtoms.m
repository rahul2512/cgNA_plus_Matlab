function Data = SugarRingAtoms(Data,pos)

% Fixed Atoms Naming and Order : O5' , O3' , C1' , N9/1 , C4/2

nbr_Seq = length(Data) ;

FixedAtomsNames = {'C1''','','','O5''','O3''','P','P',''} ;
SugarAtomsNames = {'O4''','C2''','C3''','C4''','C5'''} ;

P = diag([1 -1 -1]);

options1 = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton','MaxFunctionEvaluations',1e6,'MaxIterations',1e6);
options2 = optimset('Display','none','TolX',1e-6,'MaxFunEvals',1e4,'MaxIter',1e4);

for i = 1:nbr_Seq
    
    if nargin == 2
        
        seq = Data(i).sequence;
        
        nbp = Data(i).nbp ;
        
        bp_level = frames(Data(i).groundstate) ;
        
        [ AtomsBase, AtomsPho ] = Basepair2Atoms(bp_level, seq) ;
        
        fprintf('Reading strand progress : ')
        
        SugarCoord  = struct('W',[], 'C', cell(nbp,1)) ;
        
        SugarCoord(1).W.coord=zeros(3,6);
        SugarCoord(1).W.names= [ SugarAtomsNames {'O5'''} ];
        
        SugarCoord(1).C.coord=zeros(3,6);
        SugarCoord(1).C.names= [ SugarAtomsNames {'O5'''} ];
        
        
        for k = setdiff(2:nbp-1,pos)
            SugarCoord(k).W.coord=zeros(3,5);
            SugarCoord(k).W.names=SugarAtomsNames;
            
            SugarCoord(k).C.coord=zeros(3,5);
            SugarCoord(k).C.names=SugarAtomsNames;
        end
        
    	SugarCoord(nbp).W.coord=zeros(3,6);
        SugarCoord(nbp).W.names= [ SugarAtomsNames {'O3'''} ];
        
        SugarCoord(nbp).C.coord=zeros(3,6);
        SugarCoord(nbp).C.names= [ SugarAtomsNames {'O3'''} ];
        
        
        % =====================================================================
        %                           Watson Strand Interior
        % =====================================================================
        
        Ids = IdsPotential('int') ;
        
        j = pos;
        
        [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(j)) ;
        
        FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,j,'W');
        
        SugarCoordIni = InitializeSugarCoord(bp_level(j).Rw,FixedAtomsCoord) ;
        
        ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
        
        [ SugarCoord(j).W.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
        
        SugarCoord(j).W.names = SugarAtomsNames ;
        
        if Exit_flag == 1
            msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
        else
            [ SugarCoord(j).W.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(j).W.coord,options2) ;
            msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
        end
        

        % =====================================================================
        %                         Prepare Crick Strand data
        % =====================================================================
        
        bp_level = bp_level(nbp:-1:1);
        AtomsBase = AtomsBase(nbp:-1:1);
        AtomsPho = AtomsPho(nbp-1:-1:1);
        seq = wcc(seq,-1);
        
        % =====================================================================
        %                          Crick Strand Interior
        % =====================================================================
        
        j = pos;
        
        [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(j)) ;
        
        FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,j,'C');
        
        SugarCoordIni = InitializeSugarCoord(bp_level(j).Rc*P,FixedAtomsCoord) ;
        
        ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
        
        [ SugarCoord(j).C.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
        
        SugarCoord(j).C.names = SugarAtomsNames ;
        
        if Exit_flag == 1
            msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
        else
            [ SugarCoord(j).C.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(j).C.coord,options2) ;
            msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
        end
               
        Chain1_base = [AtomsBase.W]; Chain1_base = Chain1_base(nbp:-1:1);
        Chain1_pho  = [AtomsPho.W]; Chain1_pho = Chain1_pho(nbp-1:-1:1);
        Chain1_sugar = [SugarCoord.W];
        
        Chain2_base = [AtomsBase.C];
        Chain2_pho  = [AtomsPho.C];
        Chain2_sugar = [SugarCoord.C];
        
        Data.AtomsCoord = ArrangeAtoms( Chain1_base,Chain1_pho,Chain1_sugar, ...
            Chain2_base,Chain2_pho,Chain2_sugar  ) ;
             
    else
        
        seq = Data(i).sequence;
        
        nbp = Data(i).nbp ;
        
        bp_level = frames(Data(i).groundstate) ;
        
        [ AtomsBase, AtomsPho ] = Basepair2Atoms(bp_level, seq) ;
        
        fprintf('Reading strand progress : ')
        
        SugarCoord  = struct('W',[], 'C', cell(nbp,1)) ;
        
        % =====================================================================
        %                           Watson Strand 5' end
        % =====================================================================
        
        [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(1)) ;
        
        FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,1,'W');
        
        SugarCoordIni = InitializeSugarCoord5prime(bp_level(1).Rw,FixedAtomsCoord) ;
        
        Ids = IdsPotential('5end') ;
        
        ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
        
        [ SugarCoord(1).W.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
        
        SugarCoord(1).W.names = [ SugarAtomsNames {'O5'''} ] ;
        
        if Exit_flag == 1
            msg = [num2str(1) '(' num2str(Exit_flag) ')'] ;
        else
            
            [ SugarCoord(1).W.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(1).W.coord,options2) ;
            msg = [num2str(1) '(' num2str(Exit_flag) ')'] ;
        end
        
        fprintf( [ ' ' msg ' ' ])
        
        % =====================================================================
        %                           Watson Strand Interior
        % =====================================================================
        
        Ids = IdsPotential('int') ;
        
        for j = 2:nbp-1
            
            [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(j)) ;
            
            FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,j,'W');
            
            SugarCoordIni = InitializeSugarCoord(bp_level(j).Rw,FixedAtomsCoord) ;
            
            ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
            
            [ SugarCoord(j).W.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
            
            SugarCoord(j).W.names = SugarAtomsNames ;
            
            if Exit_flag == 1
                msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
            else
                [ SugarCoord(j).W.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(j).W.coord,options2) ;
                msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
            end
            
            fprintf( [ ' ' msg ' ' ])
            
        end
        
        % =====================================================================
        %                           Watson Strand 3' end
        % =====================================================================
        
        [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(nbp)) ;
        
        FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,nbp,'W');
        
        SugarCoordIni = InitializeSugarCoord3prime(bp_level(end).Rw,FixedAtomsCoord) ;
        
        Ids = IdsPotential('3end') ;
        
        ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
        
        [ SugarCoord(nbp).W.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
        
        SugarCoord(nbp).W.names = [ SugarAtomsNames {'O3'''} ] ;
        
        if Exit_flag == 1
            msg = [num2str(nbp) '(' num2str(Exit_flag) ')'] ;
        else
            [ SugarCoord(nbp).W.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(nbp).W.coord,options2) ;
            msg = [num2str(nbp) '(' num2str(Exit_flag) ')'] ;
        end
        
        fprintf( [ ' ' msg '\n' ])
        
        % =====================================================================
        %                         Prepare Crick Strand data
        % =====================================================================
        
        bp_level = bp_level(nbp:-1:1);
        AtomsBase = AtomsBase(nbp:-1:1);
        AtomsPho = AtomsPho(nbp-1:-1:1);
        seq = wcc(seq,-1);
        % =====================================================================
        %                           Crick Strand 5' end
        % =====================================================================
        
        [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(1)) ;
        
        FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,1,'C');
        
        SugarCoordIni = InitializeSugarCoord5prime(bp_level(1).Rc*P,FixedAtomsCoord) ;
        
        Ids = IdsPotential('5end') ;
        
        ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
        
        [ SugarCoord(1).C.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
        
        SugarCoord(1).C.names = [ SugarAtomsNames {'O5'''} ] ;
        
        fprintf('Complementary strand progress : ')
        
        if Exit_flag == 1
            msg = [num2str(1) '(' num2str(Exit_flag) ')'] ;
        else
            [ SugarCoord(1).C.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(1).C.coord,options2) ;
            msg = [num2str(1) '(' num2str(Exit_flag) ')'] ;
        end
        
        fprintf( [ ' ' msg ' ' ])
        
        % =====================================================================
        %                           Crick Strand Interior
        % =====================================================================
        
        Ids = IdsPotential('int') ;
        
        for j = 2:length(seq)-1
            
            [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(j)) ;
            
            FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,j,'C');
            
            SugarCoordIni = InitializeSugarCoord(bp_level(j).Rc*P,FixedAtomsCoord) ;
            
            ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
            
            [ SugarCoord(j).C.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
            
            SugarCoord(j).C.names = SugarAtomsNames ;
            
            if Exit_flag == 1
                msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
            else
                [ SugarCoord(j).C.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(j).C.coord,options2) ;
                msg = [num2str(j) '(' num2str(Exit_flag) ')'] ;
            end
            
            fprintf( [ ' ' msg ' ' ])
            
        end
        
        % =====================================================================
        %                           Crick Strand 3' end
        % =====================================================================
        
        [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,seq(nbp)) ;
        
        FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,nbp,'C');
        
        SugarCoordIni = InitializeSugarCoord3prime(bp_level(nbp).Rc*P,FixedAtomsCoord) ;
        
        Ids = IdsPotential('3end') ;
        
        ObjFct = @(x) Potential(x, FixedAtomsCoord, BaseId,Ids) ;
        
        [ SugarCoord(nbp).C.coord, ~, Exit_flag] = fminunc(ObjFct,SugarCoordIni,options1) ;
        
        SugarCoord(nbp).C.names = [ SugarAtomsNames {'O3'''} ] ;
        
        if Exit_flag == 1
            msg = [num2str(nbp) '(' num2str(Exit_flag) ')'] ;
        else
            [ SugarCoord(nbp).C.coord, ~, Exit_flag] = fminsearch(ObjFct,SugarCoord(nbp).C.coord,options2) ;
            msg = [num2str(nbp) '(' num2str(Exit_flag) ')' ] ;
        end
        
        fprintf( [msg  '\n'])
        
        Chain1_base = [AtomsBase.W]; Chain1_base = Chain1_base(nbp:-1:1);
        Chain1_pho  = [AtomsPho.W]; Chain1_pho = Chain1_pho(nbp-1:-1:1);
        Chain1_sugar = [SugarCoord.W];
           
        Chain2_base = [AtomsBase.C];
        Chain2_pho  = [AtomsPho.C];
        Chain2_sugar = [SugarCoord.C];
        
        Data(i).AtomsCoord = ArrangeAtoms( Chain1_base,Chain1_pho,Chain1_sugar, ...
            Chain2_base,Chain2_pho,Chain2_sugar  ) ;
        
    end
    
end


end

function [FixedAtomsNames,BaseId] = getBaseFixedAtomsName(FixedAtomsNames,base)

switch base
    
    case 'A'
        FixedAtomsNames{2} = 'N9' ;
        FixedAtomsNames{3} = 'C4' ;
        FixedAtomsNames{8} = 'C8' ;
        BaseId = 1;
        
    case 'G'
        FixedAtomsNames{2} = 'N9' ;
        FixedAtomsNames{3} = 'C4' ;
        FixedAtomsNames{8} = 'C8' ;
        BaseId = 4;
        
    case 'T'
        FixedAtomsNames{2} = 'N1' ;
        FixedAtomsNames{3} = 'C2' ;
        FixedAtomsNames{8} = 'C6' ;
        BaseId = 2;
        
    case 'C'
        FixedAtomsNames{2} = 'N1' ;
        FixedAtomsNames{3} = 'C2' ;
        FixedAtomsNames{8} = 'C6' ;
        BaseId = 3;
        
end

end

function SugarCoord = InitializeSugarCoord5prime(R_base,FixedAtomsCoord)
SugarIdealAtomsCoord;
SugarCoord = R_base*IdealSugar + FixedAtomsCoord(:,1);

O3_tmp = SugarCoord(:,8);
O3_fix = FixedAtomsCoord(:,5);
C1_fix = FixedAtomsCoord(:,1);

r1 = O3_tmp - C1_fix ;
r2 = O3_fix - C1_fix ;

theta = getValenceAngle([O3_fix C1_fix O3_tmp]);
u = cross(r1,r2);
u = u/norm(u);

Sk_u = [   0   -u(3)  u(2) ;
    u(3)   0   -u(1) ;
    -u(2)  u(1)   0 ] ;

R = cos(theta)*eye(3) + sin(theta)*Sk_u + (1-cos(theta))*(u*u');

SugarCoord = R*(SugarCoord - C1_fix) + C1_fix;

SugarCoord = SugarCoord(:,[1:5 7] );

end

function SugarCoord = InitializeSugarCoord3prime(R_base,FixedAtomsCoord)
SugarIdealAtomsCoord;
SugarCoord = R_base*IdealSugar + FixedAtomsCoord(:,1);

O5_tmp = SugarCoord(:,7) ;
C4_tmp = SugarCoord(:,4) ;
O5_fix = FixedAtomsCoord(:,4);

theta = getValenceAngle([O5_tmp C4_tmp O5_fix]);
r1 = O5_tmp - C4_tmp ;
r2 = O5_fix - C4_tmp ;

u = cross(r1,r2);
u = u/norm(u);

Sk_u = [   0   -u(3)  u(2) ;
    u(3)   0   -u(1) ;
    -u(2)  u(1)   0 ] ;

R = cos(theta)*eye(3) + sin(theta)*Sk_u + (1-cos(theta))*(u*u');

SugarCoord(:,[4 5 7]) = R*(SugarCoord(:,[4 5 6]) - C4_tmp ) + C4_tmp;

SugarCoord = SugarCoord(:,[1:5 8] );

end

function SugarCoord = InitializeSugarCoord(R_base,FixedAtomsCoord)
SugarIdealAtomsCoord;
SugarCoord = R_base*IdealSugar + FixedAtomsCoord(:,1);

O3_tmp = SugarCoord(:,8);
O3_fix = FixedAtomsCoord(:,5);
C1_fix = FixedAtomsCoord(:,1);

r1 = O3_tmp - C1_fix ;
r2 = O3_fix - C1_fix ;

theta = getValenceAngle([O3_fix C1_fix O3_tmp]);
u = cross(r1,r2);
u = u/norm(u);

Sk_u = [   0   -u(3)  u(2) ;
    u(3)   0   -u(1) ;
    -u(2)  u(1)   0 ] ;

R = cos(theta)*eye(3) + sin(theta)*Sk_u + (1-cos(theta))*(u*u');

SugarCoord = R*(SugarCoord - C1_fix) + C1_fix;

O5_tmp = SugarCoord(:,7) ;
C4_tmp = SugarCoord(:,4) ;
O5_fix = FixedAtomsCoord(:,4);

theta = getValenceAngle([O5_tmp C4_tmp O5_fix]);
r1 = O5_tmp - C4_tmp ;
r2 = O5_fix - C4_tmp ;

u = cross(r1,r2);
u = u/norm(u);

Sk_u = [   0   -u(3)  u(2) ;
    u(3)   0   -u(1) ;
    -u(2)  u(1)   0 ] ;

R = cos(theta)*eye(3) + sin(theta)*Sk_u + (1-cos(theta))*(u*u');

SugarCoord(:,[4 5 7]) = R*(SugarCoord(:,[4 5 6]) - C4_tmp ) + C4_tmp;

SugarCoord = SugarCoord(:,1:5);

end

function FixedAtomsCoord = getFixedAtomsCoord(AtomsBase,AtomsPho,FixedAtomsNames,bp_level,Strand)

nbp = length(AtomsBase);

if bp_level > 1 && bp_level < nbp
    
    names  = AtomsBase(bp_level).(Strand).names ;
    cB = AtomsBase(bp_level).(Strand).coord ;
    
    id1 = find(strncmp(names,'C1''', 3)) ;
    
    id2 = find(strncmp(names,FixedAtomsNames{2},2)) ;
    
    id3 = find(strncmp(names,FixedAtomsNames{3},2)) ;
    
    id8 = find(strncmp(names,FixedAtomsNames{8},2)) ;
    
    names = AtomsPho(bp_level-1).(Strand).names ;
    cP5 = AtomsPho(bp_level-1).(Strand).coord ;
    
    id4 = find(strncmp(names,'O5''', 3)) ;
    
    id5 = find(strncmp(names,'P', 1)) ;
    
    names = AtomsPho(bp_level).(Strand).names ;
    cP3 = AtomsPho(bp_level).(Strand).coord ;
    
    id6 = find(strncmp(names,'O3''', 3)) ;
    
    id7 = find(strncmp(names,'P', 1)) ;
    
    FixedAtomsCoord = [ cB(:,[id1 id2 id3 ]) ...
        cP5(:,id4) cP3(:,id6) ...
        cP5(:,id5) cP3(:,id7) ...
        cB(:,id8)] ; %#ok<FNDSB>
    
elseif bp_level == 1
    
    names  = AtomsBase(bp_level).(Strand).names ;
    cB = AtomsBase(bp_level).(Strand).coord ;
    
    id1 = find(strncmp(names,'C1''', 3)) ;
    
    id2 = find(strncmp(names,FixedAtomsNames{2},2)) ;
    
    id3 = find(strncmp(names,FixedAtomsNames{3},2)) ;
    
    id8 = find(strncmp(names,FixedAtomsNames{8},2)) ;
    
    names = AtomsPho(bp_level).(Strand).names ;
    cP3 = AtomsPho(bp_level).(Strand).coord ;
    
    id6 = find(strncmp(names,'O3''', 3)) ;
    
    id7 = find(strncmp(names,'P', 1)) ;
    
    FixedAtomsCoord = [ cB(:,[id1 id2 id3]) ...
        zeros(3,1) cP3(:,id6) ...
        zeros(3,1) cP3(:,id7) ...
        cB(:,id8) ] ; %#ok<FNDSB>
    
elseif bp_level == nbp
    
    
    names  = AtomsBase(bp_level).(Strand).names ;
    cB = AtomsBase(bp_level).(Strand).coord ;
    
    id1 = find(strncmp(names,'C1''', 3)) ;
    
    id2 = find(strncmp(names,FixedAtomsNames{2},2)) ;
    
    id3 = find(strncmp(names,FixedAtomsNames{3},2)) ;
    
    id8 = find(strncmp(names,FixedAtomsNames{8},2)) ;
    
    names = AtomsPho(bp_level-1).(Strand).names ;
    cP5 = AtomsPho(bp_level-1).(Strand).coord ;
    
    id4 = find(strncmp(names,'O5''', 3)) ;
    
    id5 = find(strncmp(names,'P', 1)) ;
    
    FixedAtomsCoord = [ cB(:,[id1 id2 id3]) ...
        cP5(:,id4) zeros(3,1) ...
        cP5(:,id5) zeros(3,1) ...
        cB(:,id8) ] ; %#ok<FNDSB>
    
    
end


end

function E = Potential(SugarCoord,FixedAtomsCoord,BaseId,Ids)

ForceField;

X = [ SugarCoord FixedAtomsCoord ] ;
X = X(:,Ids.X);

theta = zeros(14,1);
phi = zeros(10,1);

E = 0 ;
for i = Ids.Vangle %#ok<*NODEF>
    
    R = X(:,ValenceAngleList(i,:)) ; %#ok<*IDISVAR>
    theta(i) = getValenceAngle(R) ;
    
    p = ValenceAngleFF(ValenceAngleData(i),2:3) ;
    
    theta_0 = p{2}*pi/180;
    
    E = E + p{1}*(theta(i) - theta_0)^2 ;
    
end

for j = Ids.Dangle
    
    R = X(:,DihedralAngleList(j,:)) ;
    phi(j) = getDihedralAngle(R) ;
    
    p = DihedralAngleFF(DihedralAngleData(j),2:5) ;
    
    if size(p{1},1) > 1
        p{2} = p{2}(BaseId,:);
        p{3} = p{3}(BaseId,:);
        p{4} = p{4}(BaseId,:);
    end
    
    for k = 1:size(p{2},2)
        
        V = p{2}(k) ;
        n = p{4}(k) ;
        gamma = p{3}(k) ;
        
        E = E + V*(1+cos(n*phi(j)-gamma)) ;
        
    end
    
end

for l = Ids.NB
    
    r = X(:,NonBondedList(l,:)) ;
    r = norm(r(:,1)-r(:,2)) ;
    
    p = NonBondedFF(NonBondedData(l,:),2:3) ;
    
    R = 0.5*(p{1,1}+p{2,1}) ;
    e = sqrt(p{1,2}*p{2,2}) ;
    
    A = e*R^12; B = e*R^6;
    
    E = E + A/r^12 - B/r^6  ;
    
    q = ChargesBackboneAtoms(NonBondedList(l,:),2) ;
    
    if length(q{1})>1
        q1 = q{1}(BaseId) ;
    else
        q1 = q{1} ;
    end
    
    if length(q{2})>1
        q2 = q{2}(BaseId) ;
    else
        q2 = q{2} ;
    end
    
    E = E + q1*q2/(e*r) ;
end

for m = 1:size(BondedList,1)
    
    r = X(:,BondedList(m,:)) ;
    r = norm(r(:,1)-r(:,2)) ;
    p = BondedFF(BondedData(m),2:3) ;
    
    E = E + p{1}*(r-p{2})^2;
    
end

end

function theta = getValenceAngle(R)

r1 = R(:,1) ; r2 =R(:,2) ; r3 = R(:,3);

R1 = r1-r2;
R2 = r3-r2;

theta = acos( R1'*R2 / norm(R1) / norm(R2) )  ;

end

function phi = getDihedralAngle(R)

r1 = R(:,1) ; r2 =R(:,2) ; r3 = R(:,3); r4 = R(:,4) ;

b1 = r2 - r1 ; b2 = r3-r2 ; b3 = r4 - r3 ;

n1 = cross(b1,b2); n2 = cross(b2,b3);

m1 = cross(n1,b2/norm(b2));

x = n1'*n2 ; y = m1'*n2 ;

phi = -180/pi*atan2( y,x ) ;

end

function Ids = IdsPotential(position)

switch position
    
    case '5end'
        
        Ids.X = [1:5 7 8 9 6 11 12 13 14 ] ;
        Ids.Vangle = setdiff(1:14,11);
        Ids.Dangle = setdiff(1:15,1);
        Ids.NB = setdiff(1:9,1);
        
    case '3end'
        Ids.X = [1:5 7 8 9 10 6 12 13 14 ] ;
        Ids.Vangle = setdiff(1:14,12);
        Ids.Dangle = setdiff(1:15,[4 12]);
        Ids.NB = setdiff(1:9,7);
        
        
    case 'int'
        Ids.X = 1:13 ;
        Ids.Vangle = 1:14;
        Ids.Dangle = 1:10;
        Ids.NB = 1:9;
        
        
end

end

function AtomsCoord = ArrangeAtoms( Chain1_base,Chain1_pho,Chain1_sugar, ...
    Chain2_base,Chain2_pho,Chain2_sugar  )

nbp = length(Chain1_base);

% =========================================================================
%                               Chain 1
% =========================================================================

tmpCoord = zeros(3,5);
tmpNames = cell(1,5);

nbr_atoms = size(Chain1_base(1).coord,2);

tmpCoord(:,1:5) = Chain1_sugar(1).coord(:,[6 5 4 1 3]);
tmpNames(1:5) = Chain1_sugar(1).names(:,[6 5 4 1 3]);

tmpCoord(:,end+1) = Chain1_pho(1).coord(:,2);
tmpNames(end+1) = Chain1_pho(1).names(:,2);

tmpCoord(:,end+1) = Chain1_sugar(1).coord(:,2);
tmpNames(end+1) = Chain1_sugar(1).names(:,2);

tmpCoord(:,end+1:end+nbr_atoms) = Chain1_base(1).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain1_base(1).names(:,1:end);

AtomsCoord(1).W.coord = tmpCoord ;
AtomsCoord(1).W.names = tmpNames ;

for i=2:nbp-1
    
    tmpCoord = zeros(3,4);
    tmpNames = cell(1,4);
    
    nbr_atoms = size(Chain1_base(i).coord,2);
    
    tmpCoord(:,1:4) = Chain1_pho(i-1).coord(:,[1 4 5 3]);
    tmpNames(1:4) = Chain1_pho(i-1).names(:,[1 4 5 3]);
    
    tmpCoord(:,end+1:end+4) = Chain1_sugar(i).coord(:,[5 4 1 3]);
    tmpNames(end+1:end+4) = Chain1_sugar(i).names(:,[5 4 1 3]);
    
    tmpCoord(:,end+1) = Chain1_pho(i).coord(:,2);
    tmpNames(end+1) = Chain1_pho(i).names(:,2);
    
    tmpCoord(:,end+1) = Chain1_sugar(i).coord(:,2);
    tmpNames(end+1) = Chain1_sugar(i).names(:,2);
    
    tmpCoord(:,end+1:end+nbr_atoms) = Chain1_base(i).coord(:,1:end);
    tmpNames(end+1:end+nbr_atoms) = Chain1_base(i).names(:,1:end);
    
    AtomsCoord(i).W.coord = tmpCoord ;
    AtomsCoord(i).W.names = tmpNames ;
    
end

tmpCoord = zeros(3,4);
tmpNames = cell(1,4);

nbr_atoms = size(Chain1_base(nbp).coord,2);

tmpCoord(:,1:4) = Chain1_pho(nbp-1).coord(:,[1 4 5 3]);
tmpNames(1:4) = Chain1_pho(nbp-1).names(:,[1 4 5 3]);

tmpCoord(:,end+1:end+5) = Chain1_sugar(nbp).coord(:,[5 4 1 6 3]);
tmpNames(end+1:end+5) = Chain1_sugar(nbp).names(:,[5 4 1 6 3 ]);


tmpCoord(:,end+1) = Chain1_sugar(nbp).coord(:,2);
tmpNames(end+1) = Chain1_sugar(nbp).names(:,2);

tmpCoord(:,end+1:end+nbr_atoms) = Chain1_base(nbp).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain1_base(nbp).names(:,1:end);

AtomsCoord(nbp).W.coord = tmpCoord ;
AtomsCoord(nbp).W.names = tmpNames ;


% =========================================================================
%                               Chain 2
% =========================================================================

nbr_atoms = size(Chain2_base(1).coord,2);

tmpCoord = zeros(3,5);
tmpNames = cell(1,5);

tmpCoord(:,1:5) = Chain2_sugar(1).coord(:,[6 5 4 1 3]);
tmpNames(1:5) = Chain2_sugar(1).names(:,[6 5 4 1 3]);

tmpCoord(:,end+1) = Chain2_pho(1).coord(:,2);
tmpNames(end+1) = Chain2_pho(1).names(:,2);

tmpCoord(:,end+1) = Chain2_sugar(1).coord(:,2);
tmpNames(end+1) = Chain2_sugar(1).names(:,2);

tmpCoord(:,end+1:end+nbr_atoms) = Chain2_base(1).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain2_base(1).names(:,1:end);

AtomsCoord(1).C.coord = tmpCoord ;
AtomsCoord(1).C.names = tmpNames ;

for i=2:nbp-1
    
    tmpCoord = zeros(3,4);
    tmpNames = cell(1,4);
    
    nbr_atoms = size(Chain2_base(i).coord,2);
    
    tmpCoord(:,1:4) = Chain2_pho(i-1).coord(:,[1 4 5 3]);
    tmpNames(1:4) = Chain2_pho(i-1).names(:,[1 4 5 3]);
    
    tmpCoord(:,end+1:end+4) = Chain2_sugar(i).coord(:,[5 4 1 3]);
    tmpNames(end+1:end+4) = Chain2_sugar(i).names(:,[5 4 1 3]);
    
    tmpCoord(:,end+1) = Chain2_pho(i).coord(:,2);
    tmpNames(end+1) = Chain2_pho(i).names(:,2);
    
    tmpCoord(:,end+1) = Chain2_sugar(i).coord(:,2);
    tmpNames(end+1) = Chain2_sugar(i).names(:,2);
    
    tmpCoord(:,end+1:end+nbr_atoms) = Chain2_base(i).coord(:,1:end);
    tmpNames(end+1:end+nbr_atoms) = Chain2_base(i).names(:,1:end);
    
    AtomsCoord(i).C.coord = tmpCoord ;
    AtomsCoord(i).C.names = tmpNames ;
    
end

tmpCoord = zeros(3,4);
tmpNames = cell(1,4);

nbr_atoms = size(Chain2_base(nbp).coord,2);

tmpCoord(:,1:4) = Chain2_pho(nbp-1).coord(:,[1 4 5 3]);
tmpNames(1:4) = Chain2_pho(nbp-1).names(:,[1 4 5 3]);

tmpCoord(:,end+1:end+5) = Chain2_sugar(nbp).coord(:,[5 4 1 6 3]);
tmpNames(end+1:end+5) = Chain2_sugar(nbp).names(:,[5 4 1 6 3 ]);


tmpCoord(:,end+1) = Chain2_sugar(nbp).coord(:,2);
tmpNames(end+1) = Chain2_sugar(nbp).names(:,2);

tmpCoord(:,end+1:end+nbr_atoms) = Chain2_base(nbp).coord(:,1:end);
tmpNames(end+1:end+nbr_atoms) = Chain2_base(nbp).names(:,1:end);

AtomsCoord(nbp).C.coord = tmpCoord ;
AtomsCoord(nbp).C.names = tmpNames ;


end