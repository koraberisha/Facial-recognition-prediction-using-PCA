function [V, L, mu] = cw_pca(X)
%CW_PCA  Principal Component Analysis.
%   An example implementation for CM3106 coursework.
%
%   Input:
%     Observation matrix X. Each column of X is a datapoint (observation).
%
%   Output:
%     [V, L, mu] = kspca(X) computes the orthonormal basis V, with basis
%     vectors arranged by the corresponding eigenenergy.
%     It also returns the eigenvalues L and the mean of the data MU.


[d, N] = size(X);

% Compute the mean
mu = mean(X, 2);

% Subtract the mean from data
X = bsxfun(@minus, X, mu);

if (N >= d) % More samples than dimensions
    C = X * X'; % Covariance matrix
    [V, D, ~] = svd(C); % Singular value decomposition
else % Fewer samples than dimensions, use the trick from the remark
    C = X' * X;
    [V, D, ~] = svd(C);
    V = X * V;
    denom = repmat(sqrt(diag(D)'), d, 1);
    denom(denom == 0) = 1;
    V = V ./ denom;
end

% Eigenvalues
L = diag(D)' / N;
