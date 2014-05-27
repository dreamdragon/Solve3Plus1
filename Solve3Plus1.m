function [E,S] = Solve3Plus1(K,Kp)

% solve for the variables:
R = SolveConstrained3ptClosed(K,Kp);

R(4,:) = atan2(R(3,:),R(4,:));

S = [];
E = [];

if ~isempty(R)
    S = ones(5,size(R,2));
    S(4:5,:) = R(3:4,:);

    S(1,:) = R(2,:);
    S(2,:) = R(1,:);

    N = sqrt(sum(S(1:3,:).^2,1));
    Nm = repmat(N,3,1);
    S(1:3,:) = S(1:3,:)./Nm;
else
    return;
end;

% make essential matrices:
n = 4;
E = zeros(9,size(S,2)*n);

for i = 1:size(S,2)
    E(:,(i-1)*n+1) = reshape(skew3( S(1:3,i))*expm(skew3([0,1,0])*( S(5,i))),9,1);
    E(:,(i-1)*n+2) = reshape(skew3( S(1:3,i))*expm(skew3([0,1,0])*(-S(5,i))),9,1);
    E(:,(i-1)*n+3) = -E(:,(i-1)*n+1);
    E(:,(i-1)*n+4) = -E(:,(i-1)*n+2);
end;

