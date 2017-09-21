function [V2,Sl,Pl,Ql] = V2_Sl_calc( S,Z,V0,pm )

[~,Pt,Qt] = St(S,Z);
D = discriminant(V0,Pt,Qt);

NaN_mat = (0./(D>0)); %return nans if discriminant negative within square root

if strcmp(pm,'p') % high voltage, low losses
    V2 = ((V0^2)/2) + Pt + sqrt(D + NaN_mat);
    Slt = ((V0^2)/2) + Pt - sqrt(D + NaN_mat);
elseif strcmp(pm,'m') % low voltage, high losses
    V2 = ((V0^2)/2) + Pt - sqrt(D + NaN_mat);
    Slt = ((V0^2)/2) + Pt + sqrt(D + NaN_mat);
end

Sl = Slt./conj(Z);

Pl = real(Sl);
Ql = imag(Sl);