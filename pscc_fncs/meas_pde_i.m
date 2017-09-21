function [ I,B ] = meas_pde_i( DSSCircuit )
%MEAS_PDE_I Summary of this function goes here
%   Detailed explanation goes here

iElement = DSSCircuit.FirstPDElement;
numPDE = DSSCircuit.PDElements.Count;
i = 1;
I = nan*zeros(8,numPDE)*(1 + 1i);
B = cell(2,numPDE);
while iElement >0
    Inom = DSSCircuit.ActiveCktElement.Currents';
    Iii = Inom(1:2:end) + 1i*Inom(2:2:end);
    
    I(1:numel(Iii),i) = Iii;
    B(:,i) = DSSCircuit.ActiveCktElement.bus;
    
    iElement = DSSCircuit.NextPDElement;
    i = i+1;
end


end

