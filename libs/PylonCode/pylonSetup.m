if ~exist('QPBO.cpp','file')
    disp('QPBO-graph cut package (ver 1.3) needed to setup the pylon code.');
    disp('You can get it at Vladimir Kolmogorov''s webpage.');
    disp('Currently: http://pub.ist.ac.at/~vnk/software/QPBO-v1.3.src.tar.gz');
    disp('Please, unpack the code files into the pylon folder');
    error();
end

mex SolveQPBO.cpp QPBO_postprocessing.cpp QPBO_extra.cpp QPBO_maxflow.cpp QPBO.cpp