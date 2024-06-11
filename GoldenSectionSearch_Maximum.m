function interval = GoldenSectionSearch_Maximum (f, a, b, tol=1e-5)
  invphi = (sqrt(5) - 1) / 2;  % 1 / phi
  invphi2 = (3 - sqrt(5)) / 2; % 1 / phi^2

  [a, b] = {min(a,b), max(a,b)}{:};
  h = b - a;
  if h <= tol
    interval = [a, b];
    return; % Skip the search as tolerance is already achieved
  end

  % Required steps to achieve tolerance
  n = ceil(log(tol/h) / log(invphi));

  c = a + invphi2 * h;
  d = a + invphi * h;
  yc = f(c);
  yd = f(d);

  for k = 1:n
    if yc > yd % yc < yd to find the minimum
      b = d;
      d = c;
      yd = yc;
      h = invphi * h;
      c = a + invphi2 * h;
      yc = f(c);
    else
      a = c;
      c = d;
      yc = yd;
      h = invphi * h;
      d = a + invphi * h;
      yd = f(d);
    end
  end

  if yc < yd
    interval = [a, b];
  else
    interval = [c, d];
  end
end