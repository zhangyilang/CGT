function A = generate_A(M,N,sample_rate)
%GENERATE_A generate random Bernoulli sampling matrix

num_sample = round(sample_rate * M);
chooses = nchoosek(1:M,num_sample);

row_idxs = randperm(size(chooses,1),N);
rows = chooses(row_idxs,:);
while length(unique(rows)) < M
    row_idxs = randperm(size(chooses,1),N);
    rows = chooses(row_idxs,:);
end

idxs = sub2ind([M,N],rows,[1:N;1:N]');
A = zeros(M,N);
A(idxs) = 1;

A = A ./ sum(A,2);

end

