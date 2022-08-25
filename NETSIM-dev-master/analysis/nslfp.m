function lz = nslfp(params,v,ge,gi)

% make sure only the excitatory neurons are used
v = v( 1 : params.nbins_exc , : );
ge = ge( 1 : params.nbins_exc , : );
gi = gi( 1 : params.nbins_exc , : );

tau = 0.006 / (params.dt * params.record_downsample_factor);
exc = circshift( ge .* (params.Ee - v), tau, 2 );
inh = gi.*(params.Ei - v);
lz = zscore( abs( exc - 1.65*inh ), [], 2 );

end