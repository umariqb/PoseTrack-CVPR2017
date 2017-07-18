if (~isdeployed)
    addpath('posetrack/lib');
    addpath('posetrack/utils');
    addpath('posetrack/eval');
    addpath('deepcut/lib/pose');
    addpath('deepcut/lib/pose/multicut');
    addpath('deepcut/lib/utils');
    addpath('deepcut/lib/vis');
    addpath('deepcut/lib/multicut');
    addpath('deepcut/lib/multicut/hdf5');
    addpath('external/deepmatching');    
    caffe_dir = 'deepcut/external/caffe/matlab/';
    if exist(caffe_dir) 
        addpath(caffe_dir);
    else
        warning('Please install Caffe in ./external/caffe');
    end
    addpath('deepcut/external/liblinear-1.94/matlab/')
    fprintf('Pose startup done\n');
end
