function [y,x,A,B,C,D] = oetest2(p,u,t,x0,c)
A=[p(1),p(2),p(3),-32.17;
   p(4),p(5),p(6),p(7)
   p(8),p(9),p(10),p(11)
   0,0,1,0];
B=[p(12);p(13);p(14);p(15)];
   
[ns,ni]=size(B);
C=eye(4);
[no,ns]=size(C);
D=zeros(no,ni);
unused = c;
%
[y,x]=lsims(A,B,C,D,u,t,x0);
return