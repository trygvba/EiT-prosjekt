function bead = getBeadData(Data,column)
bead = Data(find(~isnan(Data(:,column))),column:(column+2));


end

