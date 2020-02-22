function txtldir = txtl_de_init
fp = mfilename('fullpath'); 
slashes = regexp(fp, '/');
filedir = fp(1:slashes(end)-1);
rmpath(filedir);
rmpath([filedir '/auxiliary'])
rmpath([filedir '/components'])
rmpath([filedir '/config'])
rmpath([filedir '/core'])
rmpath([filedir '/examples'])
rmpath([filedir '/examples/CompanionFiles'])
rmpath([filedir '/tests'])
rmpath([filedir '/data'])
txtldir = filedir;
end
