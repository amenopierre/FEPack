function val = FUN_per_cutoff_kite(x, vecperA, vecperB, centre, constsKite, sigma)
  % function val = FUN_per_cutoff_kite(x, vecperA, vecperB, centre, constsKite, sigma)
  % defines a 2D periodic function such that in a cell of periodicity,
  % one has a cutoff function with support inside a kite.
  %
  % The kite has the parametric equation
  %
  %       x(t) = xC + A1 * cos(t) - A2 * sin²(t)
  %       y(t) = yC + B * sin(t)
  %
  % INPUTS: * x (N x 2), the points in which the function is evaluated;
  %         * vecperA and vecperB (2 x 1), periodicity vectors;
  %         * centre (2 x 1), the centre of the cutoff
  %           centre(1) and centre(2) are reals between 0 and 1;
  %         * constKite (3 x 1), Parameters linked to the kite [A1, A2, B]
  %         * sigma plays the role of variance.

  if (nargin < 6)
    sigma = 1/sqrt(8);
  end

  if (nargin < 5)
    constsKite = [1, 1.3, 1.5]/6;
  end

  vecperA = vecperA(:);
  vecperB = vecperB(:);

  T = [vecperA, vecperB];
  invT = T \ eye(2);
  normT = norm(T, 'fro');

  Xmod = (T * (mod(invT * x(:, 1:2).', 1) - [centre(1); centre(2)])).';

  A1 = constsKite(1);
  A2 = constsKite(2);
  B  = constsKite(3);
  XmodKite = sqrt((Xmod(:, 1) + (A2 / (B^2)) * Xmod(:, 2).^2).^2 / (A1^2) +...
                  (Xmod(:, 2)).^2 / (B^2));

  val = FEPack.tools.cutoff(XmodKite, -normT, normT, sigma);
end