# Common errors running mgatk

## No samples for genotyping

#### Error message:
```
ERROR: Could not import any samples from the user specification
ERROR: check flags, logs, and input configuration (including reference mitochondrial genome)
```



## Too many open files

For certain operations, particularly in `bcall` mode, *mgatk* will attempt to have many 
file open. This can cause errors like the following:

#### Error message:
```
OSError: [Errno 24] could not open alignment file `xxxxxxxxx.bam`: Too many open files.
```

In OSX and most unix-based operating systems, the max number of file descriptors
that you can have open is fixed (for my Mac, it was 256), which can mess up the parallel 
processing performed in *mgatk*. To change this, add `ulimit -n 1024` or some
other large number to your `~/.bash_profile` as this is an environment variable that 
needs to be set in each shell session. This number should be greater than the number
of samples nominated. 

In very large data settings, it may not be possible to open as many files as there
are cells/samples (some file systems have hard upper bounds. In the `bcall` mode, 
we can get around this using the `--nsamples X` flag, where `X` is the number of
samples to be porcessed in a given setting. As long as `X` is less than the max
allowed by your OS, things should work without error. 

<br><br>
