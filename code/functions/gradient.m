function [GRAD] = gradient(ahold,y,Z,t,G,m)

    function [grad_n] = obj_n(a)

        %dimensions
        T       = size(t,2);
        N       = size(y,1)/T; 

        a=a';
        
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
            wgt(:,g) = sum(sum(((((sum(e(:,:,g).^2,2)))./((sum(e(:,:,:).^2,2)))).^(1/(m-1))),3).^(-1),3);
        end
        
        %objective function:
        
        L = sum((wgt(:,:).^m).*permute(sum(e(:,:,:).^2,2),[1,3,2]),2);
        L_n = L(n);
        
        grad_n = dlgradient(L_n,a);

    end

% matrix of gradients
for n = 1:length(y)/size(t,2)
    x0 = dlarray(ahold);
    gradval = dlfeval(@obj_n,x0);
    GRAD(:,n) =  extractdata(gradval);
end

end