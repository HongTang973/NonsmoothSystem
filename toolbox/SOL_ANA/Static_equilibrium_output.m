function equilibrium_state=Static_equilibrium_output(gap,A_1,A_2,vect_preload)
% output static points for given flight velocity
if rank(A_1)<8
    ind=1:1:6;
else
    ind=1:1:8;
end
Balance_state1=-A_1(ind,ind)\vect_preload(ind,1);
Balance_state2=-A_2(ind,ind)\vect_preload(ind,1);

if abs(Balance_state1(3))>gap
    equilibrium_state=-A_2(ind,ind)\(vect_preload(ind,1)+sign(Balance_state1(3))*vect_preload(ind,2));
%     equilibrium_state(3)=Balance_state2(3)+gap*sign(Balance_state1(3));
   
else
    equilibrium_state=Balance_state1;
end
[V1,D1] = eig(A_1);
% diag(D1)
% V1

[V2,D2] = eig(A_2);
% diag(D2)
% V2
end