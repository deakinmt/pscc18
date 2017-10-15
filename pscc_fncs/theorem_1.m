function [ Pg_prm,P0_prm ] = theorem_1( Vp, V0, Z )

lz = lambdas(Z);

Pgt_prm = Vp*(Vp - (V0*lz/sqrt(lz^2 + 1)) );
Qgt_prm = -V0*Vp/sqrt(1 + lz^2);
Sgt_prm = Pgt_prm + 1i*Qgt_prm;

Sg_prm = Sgt_prm/conj(Z);
Pg_prm = real(Sg_prm);

P0_prm = (V0^2/abs(Z))*( (Vp/V0) - (lz/sqrt(1 + lz^2)) );

end

