function val = FUN_optical_potential(x, dualVec1, dualVec2, is_even)
      
  if (is_even)
    Xmod = cos(x(:, 1:2) *  dualVec1) + ...
           cos(x(:, 1:2) *  dualVec2) + ...
           cos(x(:, 1:2) * (dualVec1  + dualVec2));

    val = FEPack.tools.cutoff(Xmod, -1.5, 1.5);
  else
    Xmod = sin(x(:, 1:2) *  dualVec1) + ...
           sin(x(:, 1:2) *  dualVec2) + ...
           sin(x(:, 1:2) * (dualVec1  + dualVec2));

    val = tanh(Xmod);
  end

end
