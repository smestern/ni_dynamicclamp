function netsim(paramfile, outputdir, jobid, options)
% Gabriel Benigno, 2021

arguments
    paramfile
    outputdir
    jobid = 0 % any nonnegative integer
    options.compile_type = 'performance' % can be 'debug', 'profiler', 'profilegen', 'profileuse', or 'performance'
    options.external_input = false 
    options.release_probability = false
    options.netsimdir = ''
end

if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end

paramfile = strrep(paramfile,' ','\ ');
outputdir = strrep(outputdir,' ','\ ');

wd = pwd;
if isempty( options.netsimdir )
    netsimdir = mfilename('fullpath');
    netsimdir = fileparts(fileparts(netsimdir));
end
cd(netsimdir)

str = sprintf('make COMPILE_TYPE=%s', options.compile_type);
if options.external_input, str = strcat(str, ' EXTERNAL_INPUT=yes'); end
if options.release_probability, str = strcat(str, ' RELEASE_PROBABILITY=yes'); end

system('make clean');
system(str);
system( sprintf('./netsim -f %s -o %s -j %u', paramfile, outputdir, jobid) ) ;

cd(wd)

end