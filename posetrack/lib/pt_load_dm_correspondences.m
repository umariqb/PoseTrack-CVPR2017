function corres_pts = pt_load_dm_correspondences(corres_fn)

try
    
if(~exist(corres_fn, 'file'))
    error('File for correspondences does not exists. %s', corres_fn);
end

corres_pts = dlmread(corres_fn);
num_cols=size(corres_pts,2);

if(num_cols < 4)
    error('Incorrect data in file: %s', corres_fn);
end

if(num_cols > 4)
    corres_pts = corres_pts(:,1:4);
end

catch 
    keyboard();
end

end