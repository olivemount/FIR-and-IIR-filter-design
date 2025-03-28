function audioobj = grabar_audio(Tc, fs, r, ID)
    audioobj = audiorecorder(fs, r, 1, ID);
    recordblocking(audioobj, Tc);
end