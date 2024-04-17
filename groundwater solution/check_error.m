function [ error_a, error_b ] = check_error( z1a,z2a,z1b,z2b,a,b )
error_a = 0;
error_b = 0;

na = 0:size(a,2)-1;
nb = 0:size(b,2)-1;

z1at = z1a;
z2at = z2a;
cnt = 0;
while ne(isempty(z1at),1)
    cnt = cnt + 1;
    z1 = z1at(1);
    pos1 = z1==z1a;
    pos2 = z1==z2a;
    pos = or(z1 == z1at,z1 == z2at);
    error_a(cnt) = sum(a(pos1,:).*cos(-na*pi),'all') - sum(a(pos2,:),'all');
    z1at(pos) = [];
    z2at(pos) = [];
end
z1bt = z1b;
z2bt = z2b;
cnt = 0;
while ne(isempty(z1bt),1)
    cnt = cnt + 1;
    z1 = z1bt(1);
    pos1 = z1==z1b;
    pos2 = z1==z2b;
    error_b(cnt) = sum(b(pos1,:).*cos(-nb*pi),'all') - sum(b(pos2,:),'all');
    pos = or(z1 == z1bt,z1 == z2bt);
    z1bt(pos) = [];
    z2bt(pos) = [];
end

error_a = abs(error_a);
error_b = abs(error_b);
end