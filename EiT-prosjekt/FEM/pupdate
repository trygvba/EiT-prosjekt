function [ p_new ] = pupdate( p,U,size )

    for i = 1:3:size
        uvec(ceil(i/3),:) = [U(i) U(i+1) U(i+2)];
    end
    p_new=p+uvec;
   
end

