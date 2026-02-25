function val = FUN_per_cutoff_cuboid(x, vecperA, vecperB, centre, BBx, BBy, sigmaX, sigmaY)
  % function val = FUN_per_cutoff_cuboid(x, vecperA, vecperB, centre, BBx, BBy, sigmaX, sigmaY)
  % defines a 2D periodic function such that in a cell of periodicity,
  % one has a cutoff function with a support in a parallelogram.
  %
  % INPUTS: * x (N x 2), the points in which the function is evaluated;
  %         * vecperA and vecperB (2 x 1), periodicity vectors;
  %         * centre (2 x 1), the centre of the cutoff
  %           centre(1) and centre(2) are reals between 0 and 1;
  %         * BBx (2 x 1), the support of the cutoff in the x-direction (centered in 0)
  %         * BBy (2 x 1), the support of the cutoff in the y-direction (centered in 0)
  %         * sigmaX plays the role of variance in the x-direction.
  %         * sigmaY plays the role of variance in the y-direction.

  if (nargin < 7)
    sigmaX = 1/sqrt(8);
  end

  if (nargin < 8)
    sigmaY = 1/sqrt(8);
  end
  
  vecperA = vecperA(:);
  vecperB = vecperB(:);

  T = [vecperA, vecperB];
  invT = T \ eye(2);
  normT = norm(T, 'fro');
  
  Xmod = (T * (mod(invT * x(:, 1:2).', 1) - [centre(1); centre(2)])).';

  val = FEPack.tools.cutoff(Xmod(:, 1), BBx(1) * normT, BBx(2) * norm(T), sigmaX) .* ...;
        FEPack.tools.cutoff(Xmod(:, 2), BBy(1) * normT, BBy(2) * norm(T), sigmaY);
        
end