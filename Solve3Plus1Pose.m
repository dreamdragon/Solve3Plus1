% Solve3Plus1Pose    Compute pose from three point correspondences K<->Kp
% and a directional correspondence g<->gp
function P = Solve3Plus1Pose(K,Kp,g,gp)

% derotate directional constraints:
c = [0,1,0]';
c = c/norm(c);

a = cross(g,c);
if norm(a) > 0.000001
    a = a/norm(a);
    theta = acos(dot(g,c));
    R = expm(skew3(a)*theta);
    Q = R*K;
    Q = Q./repmat(Q(3,:), 3,1);
else
    R = eye(3);
    Q = K;
end

ap = cross(gp,c);
if norm(ap) > 0.000001
    ap = ap/norm(ap); 
    thetap = acos(dot(gp,c));
    Rp = expm(skew3(ap)*thetap);
    Qp = Rp*Kp;
    Qp = Qp./repmat(Qp(3,:), 3,1);
else
    Rp = eye(3);
    Qp = Kp;
end

% solve:
[E,S] = Solve3Plus1(Q,Qp);

% make pose matrices:
P = zeros(16,size(S,2)*2);
sol = 1;

maxerror = 10;
for i = 1:size(S,2)
    R1 = expm(skew3(c)*( S(5,i)));
    
    t1 =   S(1:3,i);
    t2 =  -S(1:3,i);
    
    E1 = skew3(Rp'*t1)*Rp'*R1*R;
    error = max(abs(diag(Kp'*E1*K)));
    maxerror = min(maxerror,error);
    %if max(abs(diag(Kp'*E1*K))) > 0.001
    %    display('Error!!!');
        %assert(false);
    %end

    Pnew11 = [Rp'*R1*R, Rp'*t1;0,0,0,1];
    Pnew21 = [Rp'*R1*R, Rp'*t2;0,0,0,1];

    P(:,sol) = reshape(Pnew11',16,1);
    sol=sol+1;
    P(:,sol) = reshape(Pnew21',16,1);
    sol=sol+1;
end


if maxerror > 1e-6
  display('Error!!!');
        %assert(false);
end