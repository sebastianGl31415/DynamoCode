function D2=D2Matrix2ndOrder(N,h)
  D2=gallery('tridiag',N,1,-2,1)./(h^2);
end