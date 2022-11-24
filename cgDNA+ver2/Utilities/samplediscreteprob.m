function samp = samplediscreteprob(ProbMatrix)

cummat = cumsum(ProbMatrix,1);
cummat = [zeros(1,size(ProbMatrix,2)) ; cummat(1: end-1, :) ];

M = rand(1, size(ProbMatrix,2)) < cummat;
samp = sum(M,1);

end

