function proj=covproj(c,covb)
% Compute k*G*k' for multidimension G
[~,~,np]=size(covb);
[m,rho]=size(c);
if m~=np
    intermedia1=zeros(m,rho,np);
    intermedia1(:,:)=c*covb(:,:);
    intermedia2=permute(intermedia1,[2,1,3]);
    intermedia2a=intermedia2(:,:);
    intermedia3=permute(intermedia2a,[2,1])*c';
    proj=zeros(m,m,np);
    proj(:,:)=permute(intermedia3,[2,1]);
else
    % method 1
    %     proj=zeros(np,1);
    %     for ip=1:np
    %         proj(ip)=c(ip,:)*covb(:,:,ip)*c(ip,:)';
    %     end
    % method 2
    I=repmat(1:m,[rho,1]);
    trc=c';
    cmt_sp = sparse(I(:),1:m*rho,trc);
    
    I = repmat(reshape(1:rho*np,rho,1,np),[1 rho 1]);
    J = repmat(reshape(1:rho*np,1,rho,np),[rho 1 1]);
    covb_sp = sparse(I(:),J(:),covb(:));
    
    tmp_sp=cmt_sp*covb_sp;
    tmpproj=tmp_sp*cmt_sp';
    Itmp=speye(np);
    proj=tmpproj(logical(Itmp));
end
end
