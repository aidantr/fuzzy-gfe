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
    wgt(:,g) = sum(sum(((((sum(e(:,:,g).^2,2)))./((sum(e(:,:,:).^2,2)))).^(1/(m-1))),3).^(-1),3);
end

%objective function:

L = sum(sum((wgt(:,:).^m).*permute(sum(e(:,:,:).^2,2),[1,3,2]),2),1);
%L2=sum((wgt(:,:).^m).*permute(sum(e(:,:,:).^2,2),[1,3,2]),2);


%old derivative without a matrix 

%derivative wrt alpha_g
% for j=1:G
%      term1=zeros(length(y),1);
%      for g=1:G
%          if g~=j 
%              term1 = term1+((-m.*(sum((abs(del)./abs(del(:,g))).^(-2/(m-1)),2)).^(-m-1)).*((-2/(m-1)).*(((abs(del(:,j))./abs(del(:,g))).^((-m-1)/(m-1))))).*((-1.*del(:,j))./(abs(del(:,j)).*abs(del(:,g))))).*((del(:,g)).^2);
%          end
%      end
%      term2 = ((-m.*(sum((abs(del)./abs(del(:,j))).^(-2/(m-1)),2)).^(-m-1)).*(  (-2/(m-1)).*sum(   ((abs(del(:,1:G~=j)))./abs(del(:,j))).^((-m-1)/(m-1)).*((abs(del(:,1:G~=j))).*del(:,j).*((abs(del(:,j))).^(-3))),2 ))).*((del(:,g)).^2)-2.*((wgt(:,j)).^m).*del(:,j);
%      grad(j) = mean(term1+term2);
% end
% 
% %derivatives wrt time dummies and controls
% for z=1:size(Z,2)
%     term1=zeros(length(y),1);
%     term2=zeros(length(y),1);
%     for g=1:G
%         term1 = term1 + ((-m.*(sum((abs(del)./abs(del(:,g))).^(-2/(m-1)),2)).^(-m-1)).*((-2/(m-1)).*  sum( ((abs(del(:,1:G~=g)))./abs(del(:,g))).^((-m-1)/(m-1)).*(((-1.*(del(:,1:G~=g)).*Z(:,z))./((abs(del(:,1:G~=g)).*(abs(del(:,g))))))+((abs(del(:,1:G~=g)).*del(:,g).*Z(:,z))./((del(:,g)).^3))),2 ))).*((del(:,g)).^2);
%         term2 = term2 +((Wgt(:,g)).^m).*(-2).*Z(:,z).*del(:,g);
%     
%     end
%     grad(G+z) = mean(term1+term2);
% end


% 
%%%% Sandwich matrix:
%derivative wrt alpha_jtau

% for j=1:G
%     for tau=1:T
%         GRAD(:,T*(j-1)+tau) = -2.*aa.^(-m).*(sum((e(:,:,j).^2),2)).^((-m)/(m-1)).*e(:,tau,j);
%     end
% end
% 
% %derivative wrt theta_s (loop over number of common controls)
% for z=1:size(Z,2)
%     GRAD(:,G*T+z) = -2.*aa.^(-m).*sum(((sum(e(:,:,:).^2,2)).^((-m)/(m-1)).*sum(e(:,:,:).*Zt(:,:,z),2)),3);
% end

%z=1;
%test = min((sum(e(:,:,:).^2,2)).^((-m)/(m-1)),10^200);


%V = GRAD'*GRAD;

% 
%VARCOVAR = inv(H./N)*(V./N)*inv(H./N);





% create "a" vector (called aa since a is already fxn arg)
%aa      = sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1)),3);
% force not inf
%aa = max(min(aa,10^200),10^(-200));


%objective 
%L = sum(aa.^(1-m));

%g=1;
%L = sum((sum(e(:,:,g).^2,2)./(sum(e(:,:,:).^2,2))).^(1/(m-1)),3).^(1-m);

%L = sum(sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1)),3).^(1-m));
%tt = aa.^(1-m);


% %preallocate grad
% grad    = zeros(G*T+size(Z,2),1);
% 
% %derivative wrt alpha_jtau
% for j=1:G
%     for tau=1:T
%         grad(T*(j-1)+tau) = -2.*sum(aa.^(-m).*(sum((e(:,:,j).^2),2)).^((-m)/(m-1)).*e(:,tau,j));
%     end
% end
% 
% %derivative wrt theta_s (loop over number of common controls)
% for z=1:size(Z,2)
%     grad(G*T+z) = -2.* sum(aa.^(-m).*sum(((sum(e(:,:,:).^2,2)).^((-m)/(m-1)).*sum(e(:,:,:).*Zt(:,:,z),2)),3));
% end
% % 
% %preallocate hessian
% H               = zeros(G*T+size(Z,2),G*T+size(Z,2));
% 
% 
% %alpha_gt alpha_gt
% 
% for g=1:G
%     for tau=1:T
%         H(T*(g-1)+tau,T*(g-1)+tau) = 4.*sum((aa.^(-m-1)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-2*m)/(m-1))).*(e(:,tau,g).^2),1)+4.*sum((aa.^(-m)).*(-m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1) - 1)).*(e(:,tau,g).^2),1)+2.*sum((aa.^(-m)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1))),1);
%     end
% end
% 
% 
% %alpha_gt alpha_gt'
% for g=1:G
%     for tau=1:T
%         for t=1:T
%             if t~=tau
%                 H(T*(g-1)+tau,T*(g-1)+t) = 4.*sum((aa.^(-m-1)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-2*m)/(m-1))).*(e(:,tau,g).*e(:,t,g)),1)+4.*sum((aa.^(-m)).*(-m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1) - 1)).*(e(:,tau,g).*e(:,t,g)),1);
%             end
%         end
%     end
% end
% 
% %alpha_gt alpha_g't
% for g=1:G
%     for j=1:G
%         for tau=1:T
%             if g~=j
% %                H(T*(g-1)+tau,T*(j-1)+tau) = 4.*sum((aa.^(-m-1)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-2*m-2)/(m-1) - 1)).*(e(:,tau,g).*e(:,tau,j)),1);
%                 H(T*(g-1)+tau,T*(j-1)+tau) = 4.*sum((aa.^(-m-1)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-1)/(m-1) - 1)).*((sum(e(:,:,j).^2,2)).^((-1)/(m-1) - 1)).*(e(:,tau,g).*e(:,tau,j)),1);
%             end
%         end
%     end
% end
% 
% %alpha_gt alpha_g't'
% for g=1:G
%     for j=1:G
%         for tau=1:T
%             for t=1:T
%                 if g~=j & tau~=t
% %                    H(T*(g-1)+tau,T*(j-1)+t) = 4.*sum((aa.^(-m-1)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-2*m-2)/(m-1) - 1)).*(e(:,tau,g).*e(:,t,j)),1);
%                     H(T*(g-1)+tau,T*(j-1)+t) = 4.*sum((aa.^(-m-1)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-1)/(m-1) - 1)).*((sum(e(:,:,j).^2,2)).^((-1)/(m-1) - 1)).*(e(:,tau,g).*e(:,t,j)),1);
%                 end
%             end
%         end
%     end
% end
% 
% %apha_gt theta_s
% for g=1:G
%     for tau=1:T
%         for s=1:size(Z,2)
% %            H(T*(g-1)+tau,G*T+s) = (4*m/(m-1)).*sum((aa.^(-m-1)).*sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)),3).*(sum(e(:,:,g).^2,2)).^((-m)/(m-1)).*e(:,tau,g),1)-(4*m/(m-1)).*sum((aa.^(-m)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1)-1)).*sum(e(:,:,g).*Zt(:,:,s),2).*e(:,tau,g),1)+2.*sum((aa.^(-m)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1))).*Zt(:,tau,s),1);
%             H(T*(g-1)+tau,G*T+s) = (4*m/(m-1)).*sum((aa.^(-m-1)).*sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)),3).*(sum(e(:,:,g).^2,2)).^((-m)/(m-1)).*e(:,tau,g),1)-(4*m/(m-1)).*sum((aa.^(-m)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1)-1)).*sum(e(:,:,g).*Zt(:,:,s),2).*e(:,tau,g),1)+2.*sum((aa.^(-m)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1))).*Zt(:,tau,s),1);
%         end
%     end
% end
% 
% for g=1:G
%     for tau=1:T
%         for s=1:size(Z,2)
% %            H(G*T+s,T*(g-1)+tau) = (4*m/(m-1)).*sum((aa.^(-m-1)).*sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)),3).*(sum(e(:,:,g).^2,2)).^((-m)/(m-1)).*e(:,tau,g),1)-(4*m/(m-1)).*sum((aa.^(-m)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1)-1)).*sum(e(:,:,g).*Zt(:,:,s),2).*e(:,tau,g),1)+2.*sum((aa.^(-m)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1))).*Zt(:,tau,s),1);
%             H(G*T+s,T*(g-1)+tau) = (4*m/(m-1)).*sum((aa.^(-m-1)).*sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)),3).*(sum(e(:,:,g).^2,2)).^((-m)/(m-1)).*e(:,tau,g),1)-(4*m/(m-1)).*sum((aa.^(-m)).*(m/(m-1)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1)-1)).*sum(e(:,:,g).*Zt(:,:,s),2).*e(:,tau,g),1)+2.*sum((aa.^(-m)).*((sum(e(:,:,g).^2,2)).^((-m)/(m-1))).*Zt(:,tau,s),1);
%         end
%     end
% end
% 
% %theta_s theta_s'
% for s=1:size(Z,2)
%     for r=1:size(Z,2)
%         if s~=r 
% %            H(G*T+s,G*T+r) =   (4*m/(m-1)).*sum((aa.^(-m-1)).*(sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,r),2)),3)).^2,1)-(4*m/(m-1)).*sum(aa.^(-m).*(sum((sum(e(:,:,:).^2,2)).^((-m)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,r),2)).*(sum(e(:,:,:).*Zt(:,:,s),2)),3)).^2,1)+2.*sum(aa.^(-m).*(sum(sum(e(:,:,:).^2,2).^((-m)/(m-1) - 1).*(sum(Zt(:,:,s).*Zt(:,:,r),2)),3)).^2,1);
%             H(G*T+s,G*T+r) =   (4*m/(m-1)).*sum((aa.^(-m-1)).*(sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,r),2)),3)).^2,1)-(4*m/(m-1)).*sum(aa.^(-m).*(sum((sum(e(:,:,:).^2,2)).^((-m)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,r),2)).*(sum(e(:,:,:).*Zt(:,:,s),2)),3)),1)+2.*sum(aa.^(-m).*(sum(sum(e(:,:,:).^2,2).^((-m)/(m-1) - 1).*(sum(Zt(:,:,s).*Zt(:,:,r),2)),3)),1);
%         end
%     end
% end
% 
% %theta_s theta_s
% for s=1:size(Z,2)
% %       H(G*T+s,G*T+s) =  (4*m/(m-1)).*sum((aa.^(-m-1)).*(sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)),3)).^2,1)-(4*m/(m-1)).*sum(aa.^(-m).*(sum((sum(e(:,:,:).^2,2)).^((-m)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)).^2,3)).^2,1)-2.*sum(aa.^(-m).*(sum(sum(e(:,:,:).^2,2).^((-m)/(m-1) - 1).*(sum(-Zt(:,:,s).*Zt(:,:,s)+e(:,:,:),2)),3)).^2,1);
%        H(G*T+s,G*T+s) =  (4*m/(m-1)).*sum((aa.^(-m-1)).*(sum((sum(e(:,:,:).^2,2)).^((-1)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)),3)).^2,1)-(4*m/(m-1)).*sum(aa.^(-m).*(sum((sum(e(:,:,:).^2,2)).^((-m)/(m-1) - 1).*(sum(e(:,:,:).*Zt(:,:,s),2)).^2,3)),1)-2.*sum(aa.^(-m).*(sum(sum(e(:,:,:).^2,2).^((-m)/(m-1) - 1).*(sum(-Zt(:,:,s).*Zt(:,:,s)+e(:,:,:),2)),3)),1);
% end
% 


end