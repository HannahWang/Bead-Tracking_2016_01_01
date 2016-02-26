function [x,y,z,u,v,k]=trim_ty(a,b,c,d,e,f)

x = [];
y = [];
z = [];
u = [];
v = [];
k = [];

for i=1:2:13469
x=[x a(i)]; 
y=[y b(i)];
z=[z c(i)];
u=[u d(i)];
v=[v e(i)];
k=[k f(i)];

end
end