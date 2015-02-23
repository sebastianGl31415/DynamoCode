function D1=D1Matrix2ndOrder(N,h)
  D1=gallery('tridiag',N,-1,0,1)./(2*h);
end