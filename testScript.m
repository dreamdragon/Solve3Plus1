for i = 1:10000
    [Q,Qp,G,Gp,E,P3D,P,Pp] = Random3Plus1(0,5);
    if max(abs(diag(Qp'*E*Q))) > 1e-10
       display('Error!!!');
       assert(false);
    end
    P_s = Solve3Plus1Pose(Q(1:3,:),Qp(1:3,:),G,Gp);
end
