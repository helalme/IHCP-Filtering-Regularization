function [x,alf_opt] = TikhonovReg(A,b,i)

lam0=(A\b)/2;
[w,d]=eig(A);% columns of w unit norm

dist_opt=0.0;
alf_opt=0.0;
x=lam0;
g=0.00001;
%if i<=43 || i>=196
%    g=.00001;
%elseif i==44 || i==45 ||i==52 ||i==195
 %   g=1.0;
%else 
 %   g=0.1;
%end
for j=1:1
    alf=g*j;
    T(size(b))=0; %1*18
    P=T';   %18*1
    for i=1:size(b)
        temp=w(:,i);
        P = P + d(i,i)* (temp' * b * temp)/(d(i,i)*d(i,i) + alf*alf) + alf*alf*lam0/(d(i,i)*d(i,i) + alf*alf);
    end

    %selecting optimal regularization parameter
    a=log10(norm(P,2));
    b=log10(norm(A*P-b,2));
    dist_curr=sqrt((a-1)^2+(b-1)^2);
    if j==1 
        dist_opt=dist_curr;
        alf_opt=alf;
        x=P';
    end
    if j>1 
        if dist_curr<dist_opt
            dist_opt=dist_curr;
            alf_opt=alf;
            x=P';
        end
    end
end


    