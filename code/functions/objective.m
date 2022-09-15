% returns value of FCR objective function for a given vector of coefficients

function [L]=objective(a,y,Z,t,G,m)

%dimensions
T       = size(t,2);
N       = size(y,1)/T; 

%create matrix of error terms (N x T x G)
amat    = [reshape(a(1:G*T),[G,T]) repmat(a(G*T+1:end)',G,1) ];
del     = y-[t Z]*amat';
del2    = mat2cell(del,repmat(T,size(y,1)/T,1)',[G]);
e       = permute(cat(3,del2{:}),[3,1,2]);

%order controls 
Zt      = mat2cell(Z,repmat(T,size(y,1)/T,1)',[size(Z,2)]);
Zt      = permute(cat(3,Zt{:}),[3,1,2]);

%weights
wgt = zeros(N,G);

for g=1:G
    wgt(:,g) = sum(((((sum(e(:,:,g).^2,2)))./((sum(e(:,:,:).^2,2)))).^(1/(m-1))),3).^(-1);
end

%objective
L = sum(sum((wgt(:,:).^m).*permute(sum(e(:,:,:).^2,2),[1,3,2]),2),1);

end

% from the replication package for Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)
