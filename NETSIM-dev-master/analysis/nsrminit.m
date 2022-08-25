function [params,v,ge,gi,lz,st,si,ge_indiv] = nsrminit(it,params,v,ge,gi,lz,st,si,ge_indiv)
% it - initial transient duration (seconds)

DT = params.dt*params.record_downsample_factor;
n = it / DT;
m = it / params.dt;

params.t(1:n) = [];
params.t = params.t - it;
params.t = round(params.t,4);

v( : , 1:n ) = [];
ge( : , 1:n ) = [];
gi( : , 1:n ) = [];
lz( : , 1:n ) = [];
ge_indiv( : , 1:m ) = [];

cdn = (st <= it);
si(cdn) = [];
st(cdn) = [];
st = st - it;
st = round(st,4);

end