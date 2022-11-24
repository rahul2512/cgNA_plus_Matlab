load frames_2dgc.mat

tr =[  0.000000 -1.000000  0.000000      117.32000
      -1.000000  0.000000  0.000000      117.32000
       0.000000  0.000000 -1.000000       43.44000 ] ;
tr = [tr ; 0 0 0 1] ;
clc

base = base(2:end,:);
pho  = pho(2:end,:);

nbp = size(base,1);

base(end+1:2*end,:) = base ;
pho(end+1:2*end,:) = pho ; 

for i = 1:nbp
g = reshape(base(i,:) , [3 4]) ;
g = tr*[g ; 0 0 0 1] ;
g = g(1:3,:);
base(2*nbp-i+1,:) = g(:)' ;
end

for i = 1:(nbp-1)
g = reshape(pho(i,:) , [3 4]) ;
g = tr*[g ; 0 0 0 1] ;
g = g(1:3,:);
pho(2*nbp-i-1,:) = g(:)' ;
end

save('frames_2dgc_bio.mat','base','pho','tr')
