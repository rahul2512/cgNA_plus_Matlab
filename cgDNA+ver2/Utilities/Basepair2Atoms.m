function [ AtomsBase , AtomsPho ] = Basepair2Atoms( bp_level , seq )
%==========================================================================
%
% Function cgDNApdb : EmbeddedAtoms = Basepair2Atoms( basepair, Atoms )
%
% This function embed the IdealAtoms to the rigid base frames of the
% basepair structure
%
%
% Input
%      - basepair          structure basepair containing orientation and
%                          position of absolute coordinates
%      - IdealAtoms        structure containing the cartesian coordinates
%                          of the ideal atoms per each base
%
%==========================================================================



nbp = length(bp_level) ;

IdealAtoms = load('./Utilities/Ideal_Bases.mat') ;

AtomsBase  = struct('W',[], 'C', cell(nbp,1)) ;
AtomsPho  = struct('W',[], 'C', cell(nbp-1,1)) ;

j = 0; 
k = 0;

for i = 1:nbp
    base = seq(i) ;
    
    Rw = bp_level(i).Rw ; rw = bp_level(i).rw ;
    W = Rw*IdealAtoms.(base).coord + repmat(rw, [1 size(IdealAtoms.(base).coord,2) ]) ;
    
    base_c = wcc(base) ;
    Rc = bp_level(i).Rc*diag([1 -1 -1]) ; rc = bp_level(i).rc ;
    C = Rc*IdealAtoms.(base_c).coord + repmat(rc, [1 size(IdealAtoms.(base_c).coord,2) ]) ;
    
    AtomsBase(i).W.coord = W ;
    AtomsBase(i).W.names =  IdealAtoms.(base).names ;
    
    AtomsBase(i).C.coord = C ;
    AtomsBase(i).C.names =  IdealAtoms.(base_c).names ;
    
    Rpw = bp_level(i).Rpw ; rpw = bp_level(i).rpw ;
    
    if isempty(Rpw)~=1 && isempty(rpw)~=1
        
        pW = Rpw*IdealAtoms.P.coord + repmat(rpw, [1 size(IdealAtoms.P.coord,2) ]) ; 
        
        j=j+1;
  
        AtomsPho(j).W.coord = pW ;
        AtomsPho(j).W.names =  IdealAtoms.P.names ;
        
    end
    
    Rpc = bp_level(i).Rpc ; rpc = bp_level(i).rpc ;
    
    if isempty(Rpc)~=1 && isempty(rpc)~=1
        
       pC = Rpc*IdealAtoms.P.coord + repmat(rpc, [1 size(IdealAtoms.P.coord,2) ]) ; 
       
       k = k+1; 
       
       AtomsPho(k).C.coord = pC ;
       AtomsPho(k).C.names =  IdealAtoms.P.names ;
       
    end
  
    
end

end


