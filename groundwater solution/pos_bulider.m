function [ pos ] = pos_bulider( z1,z2 )

for ii = 1:length(z1)
    z1_temp = z1;
    z2_temp = z2;
    z1_is = z1(ii);
    z1_temp(ii) = [];
    z2_is = z2(ii);
    z2_temp(ii) = [];
    
    pos11 = z1_temp==z1_is;
    pos12 = z2_temp==z1_is;
    pos21 = z1_temp==z2_is;
    pos22 = z2_temp==z2_is;
    
    pos{ii} = [pos11;pos12;pos21;pos22];
end

end