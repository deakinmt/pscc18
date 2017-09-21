function Qg_p_v = qg_pv( Z,V0,Vg,Pg )

[ ~,lr,lx,~,~ ] = lambdas( Z );

lv = Vg/V0;
Sbar = (V0.^2)/abs(Z);

Qg_p_v = Sbar*( ((lv^2)*lx) - sqrt( lv^2 - ((lv.^4)*(lr.^2)) + (2*lr*(lv.^2)*(Pg/Sbar)) - ...
                                        ((Pg/Sbar).^2) )  );

                                    
[v2_idl,~,~,~] = V2_Sl_calc( Pg + 1i*Qg_p_v,Z,V0,'p' );

% to ensure that only real solutions are returned:
dV = 1e-4;
NaN_mat = 0./(sqrt(v2_idl)<Vg+dV);
Qg_p_v = Qg_p_v + NaN_mat;
                                    
end

