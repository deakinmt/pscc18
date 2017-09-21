function tot_q = calc_capQ(filename)

fileID = fopen(filename);
C = textscan(fileID,'%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
'Delimiter',',','Headerlines',1,'EmptyValue',0); %headerlines to skip first line.
fclose(fileID);

cap_q = zeros(numel(C{1}),3); %3 phases
for i = 1:numel(C{1})
    if strcmp(C{1}{i}(2:4),'Cap')
        for j=1:3
            cap_q(i,j) = C{2*j + 2}(i) + 1i*C{2*j + 3}(i);
        end
    end
end
tot_q = sum(sum(cap_q));