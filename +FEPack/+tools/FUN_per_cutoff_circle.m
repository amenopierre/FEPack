function val = FUN_per_cutoff_circle(x, vecperA, vecperB, centre, BBradius, sigma)
  % function val = FUN_per_cutoff_circle(x, vecperA, vecperB, centre, BBradius, sigma)
  % defines a 2D periodic function such that in a cell of periodicity,
  % one has a radial cutoff function.
  %
  % INPUTS: * x (N x 2), the points in which the function is evaluated;
  %         * vecperA and vecperB (2 x 1), periodicity vectors;
  %         * centre (2 x 1), the centre of the cutoff
  %           centre(1) and centre(2) are reals between 0 and 1;
  %         * BBradius (2 x 1), the support of the cutoff (centered in 0)
  %         * sigma plays the role of variance.

  if (nargin < 6)
    sigma = 1/sqrt(8);
  end

  vecperA = vecperA(:);
  vecperB = vecperB(:);

  T = [vecperA, vecperB];
  invT = T \ eye(2);
  normT = norm(T, 'fro');

  Xmod = (T * (mod(invT * x(:, 1:2).', 1) - [centre(1); centre(2)])).';
  val = FEPack.tools.cutoff(sqrt(Xmod(:, 1).^2 + Xmod(:, 2).^2), BBradius(1) * normT, BBradius(2) * normT, sigma);
end