function [psi_avg,psi_sh_avg,psi_std,psi_sh_std] = phasesimind(lz1,lz2,Fs)
% phase similarity index
% 
% lz1,lz2 -- datacubes of same size
% Fs -- sampling rate in Hz


flatten = @(x) x(:);
nt = size(lz1,3);

lz1 = bandpass_filter( lz1, 5, 50, 4, Fs );
lz2 = bandpass_filter( lz2, 5, 50, 4, Fs );

gp1 = generalized_phase( lz1, Fs, 5 );
gp2 = generalized_phase( lz2, Fs, 5 );

gp2_sh = generalized_phase( shuffle_channels(lz2) , Fs , 5 );

psi = zeros(1,nt);
psi_sh = psi;
for ii = 1:nt
    ang1 = angle(flatten(gp1(:,:,ii)));
    ang2 = angle(flatten(gp2(:,:,ii)));
    ang2_sh = angle(flatten(gp2_sh(:,:,ii)));

    psi(ii) = circ_corrcc(ang1,ang2);
    psi_sh(ii) = circ_corrcc(ang1,ang2_sh);
end

psi_avg = mean(psi);
psi_sh_avg = mean(psi_sh);

psi_std = std(psi);
psi_sh_std = std(psi_sh);

end
