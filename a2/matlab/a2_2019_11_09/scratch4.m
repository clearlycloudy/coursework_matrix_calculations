dim = 8
i=1
j=1
off=1
h=1

v = ...
-1/(h^4)*(+2*gamma(0+off,j-1+off)...
          -6*gamma(0+off,j+off)...
          +2*gamma(0+off,j+1+off))...
-1/(h^2)*(zeta(0+off,j+off))

b = ...
-1/(h^4)*(+2*gamma(i-1+off,0+off)...
          -6*gamma(i+off,0+off)...
          +2*gamma(i+1+off,0+off))...
-1/(h^2)*(eta(i+off,0+off))

e = 1/(h^4)*(2*gamma(0+off,0+off))

v+b+e

%

2*gamma(1,3)-6*gamma(1,2)+2*gamma(1,1)-6*gamma(2,1)+2*gamma(3,1)

m=1:n-1
m_minus = m-1
m_plus = m+1
m=1
u=2
arr = []
for m=2:1:n
    arr = [arr 2*gamma(1,m-1)-6*gamma(1,m)+ 2*gamma(1,m+1)];
end

arr = []
for m=2:1:n
    arr = [arr 1*gamma(1,m)];
end