function DistributionGeneration(M, dim)
% generate a certain number of gaussian distribution parameters
% @M - the total number of models, which equals to the total number of
% distributions needed
% @dim - the dimension of the multivariate distribution
gaussians = cell(M,2);
for pop = 1:M 
    % random generator
    mu = rand(1,dim);
    sigma = rand(dim,dim); 
    sigma = sigma * sigma'; % construct a positive semi-definite symmetric matrix
    gaussians{pop,1} = mu;
    gaussians{pop,2} = sigma;
end
name = '.\gaussians';
save(name, 'gaussians');

