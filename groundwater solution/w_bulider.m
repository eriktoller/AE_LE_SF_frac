function [ w ] = w_bulider( m,w,m_is,z1,z2,ab )
z1_is = z1(m_is);
z1(m_is) = [];
z2_is = z2(m_is);
z2(m_is) = [];
ab(m_is,:) = [];

n = cos((0:m-1)*pi);
nab = repmat(n,size(ab,1),1);

pos11 = z1==z1_is;
pos12 = z2==z1_is;
error_ab1 = sum(ab(pos11,:).*nab(pos11,:),'all') - sum(ab(pos12,:),'all');
pos21 = z1==z2_is;
pos22 = z2==z2_is;
error_ab2 = sum(ab(pos21,:).*nab(pos21,:),'all') - sum(ab(pos22,:),'all');

w(end+1) = -error_ab1;
w(end+1) = error_ab2;

end