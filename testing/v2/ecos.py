A
n = 10; m = 3; density = 0.50;  
  
% linear term
c = [zeros(n,1); ones(n,1)];  
  
% equality constraints
A = sprandn(m,n,density);  
Atilde = [A, zeros(m,n)];  
b = randn(m,1);  
  
% linear inequality constraints
I = speye(n);  
G = [  I -I;  
      -I -I];  
h = zeros(2*n,1);  
  
% cone dimensions (LP cone only)
dims.l = 2*n; dims.q = [];  
  
% call solver
z = ecos(c,G,h,dims,Atilde,b);  
x = z(1:n);  
nnzx = sum(abs(x) > 1e-8);  
  
% print sparsity info
fprintf('Optimal x has %d/%d (%4.2f%%) non-zero entries.\n', nnzx , n,  nnzx/n*100);  




# apparently you can feed it a first and second order cone.
# so the first order cone is for linear constraints
# the second order cone is for the quadratic somehow
