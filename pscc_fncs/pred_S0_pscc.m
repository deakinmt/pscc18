function [ S0,Sg,Sl,V2,I,P0,Q0 ] = pred_S0_pscc( Pgen, Sload, Z, V0, Vp  )

Pg = Pgen - real(Sload);
Qg = real(qg_pv( Z,V0,Vp,Pg )); %avoid any complex artefacts
Sg = Pg + 1i*Qg;

[ V2,Sl ] = V2_Sl_calc( Sg,Z,V0,'p' );

P0 = Pg - real(Sl);
Q0 = Qg - imag(Sl);
S0 = P0 + 1i*Q0;

I = sqrt(abs(Sl/Z)); %pu

end

