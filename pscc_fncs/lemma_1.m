function [ Pghat, P0hat ] = lemma_1( Vp, Ip, V0, Z )

Pgt_hat = 0.5*(Vp^2 - V0^2 + abs((Ip*Z))^2);
Qgt_hat = -sqrt(abs(Vp*Ip*Z)^2 - Pgt_hat^2);

Sghat = (Pgt_hat + 1i*Qgt_hat)/conj(Z);
Pghat = real(Sghat);

P0hat = Pghat - real(Z)*(Ip^2);

end