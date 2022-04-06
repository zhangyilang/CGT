function [x,xbits,x_pos] = GenSparseVec_COVID19(lambda,mu,sigma,N)
%GENSPARSEVEC_COVID19 generate sparse vectors from uniform distribution
% the occurance of non-zero xi ~ Poisson(lambda)
% each non-zero element ln xi ~ N(mu,sigma)
% N     length of the vector
% K     sparsity of the vector

x               = zeros(N,1);
xbits           = zeros(N,1);

K               = poissrnd(lambda * N);
x_pos           = randperm(N,K);

x(x_pos)        = 10.^(normrnd(mu,sigma,[K,1]));
xbits(x_pos)    = 1;

end

