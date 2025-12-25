function parsave(fname, x, is_struct, varargin)

  if (is_struct)
    save(fname, '-struct', 'x', '-v7.3');
  else
    save(fname, 'x', '-v7.3');
  end
end