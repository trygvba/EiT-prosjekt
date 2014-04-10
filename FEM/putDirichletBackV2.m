function U = putDirichletBackV2(Uprev, lowerPlate, upperPlate, uzilow, uziup, x_plates, y_plates)

data = [3*upperPlate uziup; 3*lowerPlate uzilow; (3*(x_plates-1)+1) zeros(length(x_plates),1); (3*(y_plates-1)+2) zeros(length(y_plates),1)];
nodes = sort([3*upperPlate; 3*lowerPlate;(3*(x_plates-1)+1); (3*(y_plates-1)+2)],'ascend');
U = Uprev;
n=length(nodes);
for j=1:n
    i = nodes(j);
    temp = data(find(data(:,1)==i),2);
    N = length(U);
    if i==1
        U = [temp; U];
    elseif i==N
        U = [U; temp];
    else
        U = [U(1:(i-1)); temp; U(i:N)];
    end
end


end