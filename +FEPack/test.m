% Example matrix
n = 10;
B = diag(1:n);    % eigenvalues 1,...,n
Afun = @(lambda) (lambda.^2)*eye(n) - B;

% Scan condest over a grid
lam = linspace(-5,5,400);
cvals = zeros(size(lam));
for k = 1:numel(lam)
    try
        cvals(k) = condest(Afun(lam(k)));
    catch
        cvals(k) = Inf;   % singular matrix
    end
end

% Plot
semilogy(lam, cvals), grid on
xlabel('\lambda'), ylabel('condest(A(\lambda))')
title('Blow-ups indicate nonlinear eigenvalues')