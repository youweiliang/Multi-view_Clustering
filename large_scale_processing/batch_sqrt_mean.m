function mu = batch_sqrt_mean(W, n_chunk)
batch_size = ceil(length(W) / n_chunk);
batch_mu = cell(1, n_chunk);
for i=1:n_chunk
    first_idx = (i-1)*batch_size + 1;
    last_idx = min(i*batch_size, length(W));
    batch_mu{i} = mean(sqrt(W(first_idx:last_idx, :)), 2);
end
mu = cat(1, batch_mu{:});
mu = mean(mu);