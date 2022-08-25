function makeParameterFile(paramfile, options)

arguments
    paramfile
    
    options.N {mustBeInteger,mustBePositive}
    options.K {mustBeInteger,mustBePositive}
    options.T {mustBeNumeric,mustBePositive}
    options.ge {mustBeNumeric,mustBePositive}
    options.gi {mustBeNumeric,mustBePositive}
    options.connector {mustBeInteger}
    options.L {mustBeNumeric,mustBePositive}
    options.sigma_space {mustBeNumeric,mustBePositive}
    
    options.initiator {mustBeInteger} = 0
    options.min_conduction_speed {mustBeNumeric,mustBePositive} = 0.2
    options.max_conduction_speed {mustBeNumeric,mustBePositive} = 0.2
    options.dt {mustBeNumeric,mustBePositive} = 1e-4
    options.taum {mustBeNumeric,mustBePositive} = 20e-3
    options.vr {mustBeNumeric} = -70e-3
    options.vreset {mustBeNumeric} = -70e-3
    options.vth {mustBeNumeric} = -50e-3
    options.taur {mustBeNumeric,mustBePositive} = 5e-3
    options.Ee {mustBeNumeric} = 0
    options.Ei {mustBeNumeric} = -80e-3
    options.taue {mustBeNumeric,mustBePositive} = 5e-3
    options.taui {mustBeNumeric,mustBePositive} = 5e-3
    options.El {mustBeNumeric} = -65e-3
    options.Ie {mustBeNumeric} = 0
    options.Cm {mustBeNumeric,mustBePositive} = 200e-12
    options.synapse_delay {mustBeNumeric,mustBePositive} = 3e-4
    options.start_record_time {mustBeNumeric,mustBeNonnegative} = 0
    options.stop_record_time {mustBeNumeric,mustBePositive}
    options.bin_size {mustBeInteger,mustBePositive} = 10
    options.record_downsample_factor {mustBeInteger,mustBePositive} = 2
    options.save_connectivity {mustBeInteger} = 0
    options.vm_mean {mustBeNumeric} = -0.0680901
    options.vm_sigma {mustBeNumeric} = 0.0057296
    options.ge_mean {mustBeNumeric} = 1.14244e-08
    options.ge_sigma {mustBeNumeric} = 4.30869e-09
    options.gi_mean {mustBeNumeric} = 7.95233e-08
    options.gi_sigma {mustBeNumeric} = 5.80938e-08
    options.report_minutes {mustBeNumeric,mustBePositive} = .5
    options.p_release = 1
end

if ~exist( fileparts(paramfile) , 'dir' )
    mkdir(fileparts(paramfile))
end

if ~isfield(options,'stop_record_time')
    options.stop_record_time = options.T;
end

% read template
wd = fileparts( mfilename('fullpath') );
templ = sprintf('%s/templates/template.parameters',wd);
fid = fopen(templ, 'r');
f = fread(fid, '*char')';
fclose(fid);

% integers
f = parseWriteInts(f,options.N,'N');
f = parseWriteInts(f,options.K,'K');
f = parseWriteInts(f,options.connector,'connector');
f = parseWriteInts(f,options.initiator,'initiator');
f = parseWriteInts(f,options.bin_size,'bin_size');
f = parseWriteInts(f,options.record_downsample_factor,'record_downsample_factor');
f = parseWriteInts(f,options.save_connectivity,'save_connectivity');

% floats
f = parseWriteFloats(f,options.T,'T');
f = parseWriteFloats(f,options.ge,'ge');
f = parseWriteFloats(f,options.gi,'gi');
f = parseWriteFloats(f,options.L,'L');
f = parseWriteFloats(f,options.sigma_space,'sigma_space');
f = parseWriteFloats(f,options.min_conduction_speed,'min_speed');
f = parseWriteFloats(f,options.max_conduction_speed,'max_speed');
f = parseWriteFloats(f,options.dt,'dt');
f = parseWriteFloats(f,options.taum,'taum');
f = parseWriteFloats(f,options.vr,'vr');
f = parseWriteFloats(f,options.vreset,'vreset');
f = parseWriteFloats(f,options.vth,'vth');
f = parseWriteFloats(f,options.taur,'taur');
f = parseWriteFloats(f,options.Ee,'Ee');
f = parseWriteFloats(f,options.Ei,'Ei');
f = parseWriteFloats(f,options.taue,'taue');
f = parseWriteFloats(f,options.taui,'taui');
f = parseWriteFloats(f,options.El,'El');
f = parseWriteFloats(f,options.Ie,'Ie');
f = parseWriteFloats(f,options.Cm,'Cm');
f = parseWriteFloats(f,options.synapse_delay,'synapse_delay');
f = parseWriteFloats(f,options.start_record_time,'start_record_time');
f = parseWriteFloats(f,options.stop_record_time,'stop_record_time');
f = parseWriteFloats(f,options.record_downsample_factor,'record_downsample_factor');
f = parseWriteFloats(f,options.vm_mean,'vm_mean');
f = parseWriteFloats(f,options.vm_sigma,'vm_sigma');
f = parseWriteFloats(f,options.ge_mean,'ge_mean');
f = parseWriteFloats(f,options.ge_sigma,'ge_sigma');
f = parseWriteFloats(f,options.gi_mean,'gi_mean');
f = parseWriteFloats(f,options.gi_sigma,'gi_sigma');
f = parseWriteFloats(f,options.report_minutes,'report_minutes');
f = parseWriteFloats(f,options.p_release,'p_release');

% write output
fid = fopen(paramfile, 'w');
fprintf(fid, '%s', f);
fclose(fid);

end




function f = parseWriteInts(f,value,var)
val = sprintf( '%d', value );
expr = sprintf('<<%s>>',var);
f = strrep(f, expr, val);
end

function f = parseWriteFloats(f,value,var)
val = sprintf( '%.10e', value );
expr = sprintf('<<%s>>',var);
f = strrep(f, expr, val);
end
