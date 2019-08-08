function[Xh,U,Y,gamma]=mypca(X,p)
d=size(X,1); n=size(X,2);
X=X-repmat(mean(X),d,1);

R=X*X';
[U,gamma]=eigs(R,p);
Y=U'*X;
Xh=U*Y;

end
 