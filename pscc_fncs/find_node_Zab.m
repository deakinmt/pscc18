function Zab = find_node_Zab( Zbus,idxA,idxB )

Zab=zeros(numel(idxA),numel(idxB));

for i = 1:numel(idxA)
    for j = 1:numel(idxB)
            Zab(i,j) = Zbus(idxA(i),idxB(j));
    end
end

end