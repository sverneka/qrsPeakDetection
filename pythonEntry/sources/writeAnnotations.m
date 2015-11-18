function writeAnnotations(file)
    peaks = load('../peaks.txt');
    [isloaded,config] = wfdbloadlib;
    if(config.inOctave)
    crash_dumps_octave_core(0);
    end
    fprintf('isloaded = %d\n',isloaded);
    annName = 'qrs';
    wrann(file,annName,peaks,[],[],[],[]);
end
