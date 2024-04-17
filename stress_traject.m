function [z_lines1,z_lines2] = stress_traject(space,dim,ds,Num,lvs,sigma_func)
% Fucntion settings
cond = 1e-5; % codition number for numerical iterator to exit
scale = 1; % amount to scale the levels for even value spacing, 1=equal space

% Get the starting polints
if space == 1
    % Evenly spaced starting potins
    zz = linspace(dim(1),dim(2),lvs);
    zz2 = linspace(dim(3),dim(4),lvs);
else
    % Evenly spaced values for starting potins
    lvs_calc = lvs*scale+1;
    zzt = linspace(dim(1),dim(2),lvs_calc);
    zz2t = linspace(dim(3),dim(4),lvs_calc);
    s1 = zeros(1,lvs_calc); s2 = s1;
    for ii = 1:length(zzt)
        s1(ii) = sigma_func(zzt(ii));
        s2(ii) = sigma_func(zz2t(ii));
    end
    del_s1 = max(s1)-min(s1);
    del_s2 = max(s2)-min(s2);
    zz = []; zz2 = zz;
    for ii = 1:lvs
        diffs1 = abs(s1 - max(s1) + del_s1/lvs*(ii-1));
        pos1 = find(diffs1-min(diffs1)<cond);
        zz = [zz,zzt(pos1)];
        diffs2 = abs(s2 - max(s2) + del_s2/lvs*(ii-1));
        pos2 = find(diffs2-min(diffs2)<cond);
        zz2 = [zz2,zz2t(pos2)];
    end
end

%% TEMP
zz2t = complex(zz,dim(3));
zzt = complex(dim(1),zz2);
zz = zzt;
zz2 = zz2t;

% Estimate the time
tic
for ii = 1:1
    % Get the lines for the odd numbers, i.e. goign from - to +
    % Get the loop value
    kk = ii;
    % Get the simga_1 lines
    z = zz(kk);
    zold = z;
    ds1 = ds;
    for jj = 2:Num
        [sigma_1,~,theta_p] = sigma_func(z);
        sigma = abs(sigma_1)*exp(1i*theta_p);
        z11 = z + conj(sigma)/abs(sigma)*ds1;
        if norm(z11-zold) < sqrt(2*ds^2) && jj > 4
            ds1 = -ds1;
        end
        zrk = z + conj(sigma)/abs(sigma)*ds1;
        
        ee = 1;
        NIT = 0;
        while ee > cond && NIT < 10
            zrk_old = zrk;
            [sigma_1,~,theta_p] = sigma_func(zrk);
            sigma1 = abs(sigma_1)*exp(1i*theta_p);
            zrk = z + conj(sigma + sigma1)/abs(sigma + sigma1)*ds1;
            ee = norm(zrk_old - zrk);
            NIT = NIT + 1;
        end
        
        zold = z;
        z = zrk;
    end

    % Get the simga_1 lines
    z = zz2(kk);
    zold = z;
    ds1 = ds;
    for jj = 2:Num
        [~,sigma_2,theta_p] = sigma_func(z);
        sigma = abs(sigma_2)*exp(1i*(theta_p+pi/2));
        z11 = z + conj(sigma)/abs(sigma)*ds1;
        if norm(z11-zold) < sqrt(2*ds^2) && jj > 4
            ds1 = -ds1;
        end
        zrk = z + conj(sigma)/abs(sigma)*ds1;
        
        ee = 1;
        NIT = 100;
        while ee > cond && NIT < 10
            zrk_old = zrk;
            [~,sigma_2,theta_p] = sigma_func(zrk);
            sigma1 = abs(sigma_2)*exp(1i*(theta_p+pi/2));
            zrk = z + conj(sigma + sigma1)/abs(sigma + sigma1)*ds1;
            ee = norm(zrk_old - zrk);
            NIT = NIT + 1;
        end
        

        zold = z;
        z = zrk;
    end
end
% Creating the progress bar and estimating the time left
time = toc;
formatSpec = 'Time remining: %d sec';
str = sprintf(formatSpec,ceil(time*(lvs*4-ii)));
f = waitbar(0,str,'Name','Plotting contour');

% Produces the empty plot matricies
lvs_z1 = length(zz);
lvs_z2 = length(zz2);
z_lines1 = zeros(lvs_z1,Num);
z_lines2 = zeros(lvs_z2,Num);

for ii = 1:lvs_z1
    % Get the simga_1 lines going +ds start
    z = zz(ii);
    z_lines1(ii,1) = z;
    zold = z;
    ds1 = ds;
    for jj = 2:Num
        [sigma_1,~,theta_p] = sigma_func(z);
        sigma = abs(sigma_1)*exp(1i*theta_p);
        z11 = z + (sigma)/abs(sigma)*ds1; %no conj
        zrk = z + (sigma)/abs(sigma)*ds1;%no conj
        
        ee = 1;
        NIT = 0;
        while ee > cond && NIT < 10
            zrk_old = zrk;
            [sigma_1,~,theta_p] = sigma_func(zrk);
            sigma1 = abs(sigma_1)*exp(1i*theta_p);
            zrk = z + (sigma + sigma1)/abs(sigma + sigma1)*ds1;%no conj
            ee = norm(zrk_old - zrk);
            NIT = NIT + 1;
        end
        
        z_lines1(ii,jj) = zrk;
        zold = z;
        z = zrk;
    end

    str = sprintf(formatSpec,ceil(time*(lvs_z1+lvs_z2-ii)));
    waitbar(ii/(lvs_z1+lvs_z2),f,str)
end

for ii = 1:lvs_z2
    % Get the simga_2 lines going +ds start
    z = zz2(ii);
    z_lines2(ii,1) = z;
    zold = z;
    ds1 = ds;
    for jj = 2:Num
        [~,sigma_2,theta_p] = sigma_func(z);
        sigma = abs(sigma_2)*exp(1i*(theta_p+pi/2));
        zrk = z + (sigma)/abs(sigma)*ds1; %no conj
        
        ee = 1;
        NIT = 0;
        while ee > cond && NIT < 10
            zrk_old = zrk;
            [~,sigma_2,theta_p] = sigma_func(zrk);
            sigma1 = abs(sigma_2)*exp(1i*(theta_p+pi/2));
            zrk = z + (sigma + sigma1)/abs(sigma + sigma1)*ds1; %no conj
            ee = norm(zrk_old - zrk);
            NIT = NIT + 1;
        end
        z_lines2(ii,jj) = zrk;
        zold = z;
        z = zrk;
    end
    
    str = sprintf(formatSpec,ceil(time*(lvs_z2-ii)));
    waitbar((ii+lvs_z1)/(lvs_z1+lvs_z2),f,str)
end

close(f)

end

