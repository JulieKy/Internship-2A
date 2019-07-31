[x Fs] = audioread('1.mp3');
    %% resampling to 4000 Hz
    xs=resample(x,4000,Fs);
    fn=4000;
    %% filtering BP 100-1000Hz
    y = filterbp(xs,fn);
[output_temporal_features] = temporal_features(y,fn);