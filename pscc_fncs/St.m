function [ST,PT,QT] = St(S,Z)

ST = S.*conj(Z);
PT = real(ST);
QT = imag(ST);

end

