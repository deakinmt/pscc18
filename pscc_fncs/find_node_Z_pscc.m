function [ ZZ ] = find_node_Z_pscc( YZNodeOrder,Zbus,bus1,bus2,nTr )

idx1 = find_node_idx(YZNodeOrder,bus1);
idx2 = find_node_idx(YZNodeOrder,bus2);

Z11 = find_node_Zab( Zbus,idx1,idx1 );
Z12 = find_node_Zab( Zbus,idx1,idx2 );
Z21 = find_node_Zab( Zbus,idx2,idx1 );
Z22 = find_node_Zab( Zbus,idx2,idx2 );

aa = exp(1i*2*pi/3);

MM = -nTr*sqrt(1/3)*[1; aa^2; aa]; % scaled to make it orthonormal
NN = +sqrt(1/3)*[1; aa^2; aa];

dV = NN'*( (Z21*MM) + (Z22*NN) ) - (-MM)'*( (Z11*MM) + (Z12*NN)) ;
dI = 1;
ZZ = dV/dI;

end

