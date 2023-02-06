# lbsim

To build:
```
mkdir build
cd build
cmake ..
make
```
## Input

To use with simulated load:
```
./lbsim <num PEs: default 8192> <num objs: default 65536> <RNG seed: default -1> -json <path to JSON load config file>
```

To use with VT logs (note that a _phase_ in VT parlance is what we'd call an _iteration_ in Charm++ parlance):
```
./lbsim -phase <phase to use from VT logs> -vt <paths to VT load log files..., e.g. vt_log_files/*>
```

## Output
Prints a block for each LB like so:

```
Elapsed time for rkd4: 0.478792
Maxloads: 172.403 246.823 (∑=419.225, max=246.823)
Ratio: 1.05901 1.01116 (∑=1.03031, max=1.01116)
```
`Maxloads` gives the maximum load in each dimension for the mapping, parenthetically followed by the sum and maximum across the dimensions.

`Ratio` gives the max/avg ratio in each dimension for the mapping, parenthetically followed by the ratios for the sum and maximum scoring functions.
(N.B.: The ratio for the maximum scoring function is calculated via (max over all PE load vectors) / (∑ of maximum value in load vector of each PE)).
