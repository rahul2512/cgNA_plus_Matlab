pos = 13;
dist =7;
nsnap = 10000;

res = nan(2*dist,nsnap);

%seq = 'GCTTAGTTCAAATTTGAACTAAGC' ;
seq = 'GCTCTCTGTATTAATACAGAGAGC';
%seq = repmat('AT',[1 12]);
%seq = repmat('G',[1 24]);

Data = cgDNAp(seq);

C = inv(Data.stiff);
C = (C+C')*0.5;
mu = Data.groundstate;
w = mvnrnd(mu,C,nsnap)';


for i = 1:nsnap
   
   bp_level = frames(w(:,i)); 
   
   Pw = bp_level(pos).rpw;
   
   for j = 1:dist
   
       Pc_up   = bp_level(pos+j-1).rpc;
       Pc_down = bp_level(pos-j).rpc;
   
       res(dist+j,i) = norm(Pc_up-Pw);
       res(dist-j+1,i) = norm(Pc_down-Pw);
       
   end
       
   minor(i) = min(res(dist:end,i));
   major(i) = min(res(1:dist-1,i));
   
end

figure 
imagesc(res)
colormap('Bone')
% histogram(minor);
% hold on
% histogram(major);

hold off 