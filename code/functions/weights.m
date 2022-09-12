function weight = weights(alpha,zeta,y,t,Z,G,m)

T=size(t,2);
N = size(y,1)/T; 

amat=[reshape(alpha,[G,T]) repmat(zeta',G,1)];
%amat = [reshape(fe,[T,group])' repmat(td',group,1) repmat(zeta',group,1)];
del  = y-[t Z]*amat';
del2 = mat2cell(del,repmat(T,size(y,1)/T,1)',[G]);
e = permute(cat(3,del2{:}),[3,1,2]);

%weights
wgt = zeros(N,G);
for g=1:G
    wgt(:,g) = sum(((((sum(e(:,:,g).^2,2))./((sum(e(:,:,:).^2,2)))).^(1/(m-1)))),3).^(-1);
end

weight = repelem(wgt,T,1);

end

% from the replication package for Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)