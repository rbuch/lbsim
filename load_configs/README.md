The load configuration file specifies the computational load, the data size of the chare array, and the size of messages sent at the end of phases.

There are three top level keys:

1. "phases" - An array that specifies load distributions for each phase.
2. "dataSize" - An array of size one that specifies the distribution of the data size of each object.
3. "commSizes" - An array of the same size as the "phases" array that specifies the size of the messages sent to neighbors after each phase. Currently this only supports using ints directly, not distributions.

There are several available distributions:

1. "constant" - A constant distribution of a single value.
       - "value" - The value to use.
2. "linear" - A linear distribution of the form `base` + `increment` * `index`.
   `index` is calculated by taking the 3D index of the chare array object and
   flattening it. `shift` allows the index of the base to be modularly shifted
   by the specified value (e.g. a `shift` of 1 would change {`base`, `base` +
   `increment`, `base` + `increment` * 2} to {`base` + `increment` * 2, `base``,
   `base` + `increment`}).
      - "base" - The base value.
      - "increment" - The value to increment by.
      - "shift" - Number of elements to shift the distribution by. Defaults to 0.
3. "normal" - A normal distribution. Note that negative loads are not allowed, so any negative samples are replaced with 0.
      - "mean" - The mean.
      - "stddev" - The standard deviation.
3. "exponential" - An exponential distribution.
      - "lambda" - The lambda parameter.
3. "nestedProb" - A container for probabilistically selecting one of several distributions.
      - "ratio" - An array of ints specifying the ratios for selecting a distribution. `[1, 1]` selects between two distributions, each with probablility 0.5, while `[2, 1, 1]` selects from three distributions with probablility 0.5, 0.25, 0.25.   
      - "dists" - An array of distributions to select from, which can be any of available distributions. The size of this array must match the size of the "ratio" array. 
4. "nestedBlock" - A container for selecting one of several distributions in a blocked manner.
      - "ratio" - An array of ints specifying the ratios for selecting a distribution. `[1, 1]` assigns the first distribution to the first half of the elements and the second distribution to the second half, while `[2, 1, 1]` assigns the first distribution to the first half of the elements, the second distribution to the third quarter of the elements, and the third distribution to the last quarter of the elements.   
      - "dists" - An array of distributions to select from, which can be any of available distributions. The size of this array must match the size of the "ratio" array. 

