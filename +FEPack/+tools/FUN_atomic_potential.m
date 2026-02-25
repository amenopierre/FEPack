function val = FUN_atomic_potential(x, vecPer1, vecPer2, centers, amps, rads)

  % centers: 2 x Nc vector
  % amps:   Nc x 1  vector
  % rads:   Nc x 1  vector
  if (nargin < 6), rads = 0.4; end
  if (nargin < 5), amps = 10; end
  if (nargin < 4), centers = [0; 0]; end

  vecPer1 = vecPer1(:);
  vecPer2 = vecPer2(:);

  R = [vecPer1, vecPer2];
  T = R \ eye(2);
  normR = norm(R, 'fro');

  Nc = size(centers, 2);
  val = zeros(size(x, 1), 1);

  for idC = 1:Nc

    Xmod = (R * (mod(T * (x(:, 1:2).' - centers(:, idC)) + 0.5, 1) - 0.5)).';
    normXmod = sqrt(Xmod(:, 1).^2 + Xmod(:, 2).^2);

    val = val + amps(idC) * FEPack.tools.cutoff(normXmod, -rads(idC)*normR, rads(idC)*normR);

  end

end