# genotype

[![NPM](https://nodei.co/npm/genotype.png?global=true)](https://nodei.co/npm/genotype/)

functions for manipulating genotypes and genotype likelihoods

## Overview

Genotype likelihoods describe the probability of a set of data (typically sequencing observations) given a particular genotype. They are the common currency of genotype-based inference. This module encodes a number of likelihood methods in order to enable their comparison in a common, simple framework.

The module exports a Genotype and Observation class. Genotypes are encoded as multisets of alleles, where the alleles can be objects of any hashable type. Observations link an allele, a count, and a quality, which is assumed to be in [phred](https://en.wikipedia.org/wiki/Phred_quality_score) format.

## Installation

```
npm install genotype
```

## Usage

``` js
var Genotype = require('genotype').Genotype
console.log(new Genotype([0,1]))
/*
{ alleles: [ 0, 1 ],
  ploidy: 2,
  alleleCount: { '0': 1, '1': 1 },
  homozygous: false }
*/
var Observation = require('genotype').Observation
console.log(new Observation(1, 34))
// { allele: 1, quality: -7.828789316179756 }
```
Input qualities are converted to log10, which is used internally for calculations.

## Testing genotype likelihoods

The executable `cli.js` provides a facility to test the genotype likelihood models from [freebayes](https://github.com/ekg/freebayes) and [samtools](http://www.htslib.org/). A hard-coded quality distribution is encoded in the module, but you can vary the depth of coverage maximum and the number of trials. For each trial, a depth is sampled uniformly from `[1,depth]` and reads are generated from an underlying (hidden) genotype that is uniformly sampled from 0|0, 0|1, 1|0, 1|1. Errors are simulated by sampling an observation quality from a predefined quality distribution and inducing an error with probability proportional to the error probability represented quality score. The output is tsv-formatted and has a header which describes each field.

``` bash
% node cli.js 10000 20 >results.tsv
% head results.tsv
method  true.genotype   status  called.genotype true.genotype.likelihood        n.observations  p.homref        p.het   p.homalt        observations
freebayes       0/0     pass    0/0     0       8       0       -7.059117788843288      -24.177143476437486     [{"allele":1,"quality":{"quality":4},"error":true},{"allele":0,"quality":{"quality":2},"error":false},{"allele":0,"quality":{"quality":12},"error":false},{"allele":0,"quality":{"quality":23},"error":false},{"allele":0,"quality":{"quality":12},"error":false},{"allele":0,"quality":{"quality":17},"error":false},{"allele":0,"quality":{"quality":18},"error":false},{"allele":0,"quality":{"quality":25},"error":false}]
samtools        0/0     pass    0/0     0       8       0       -1.4872059281142307     -24.177143476437486     [{"allele":1,"quality":{"quality":4},"error":true},{"allele":0,"quality":{"quality":2},"error":false},{"allele":0,"quality":{"quality":12},"error":false},{"allele":0,"quality":{"quality":23},"error":false},{"allele":0,"quality":{"quality":12},"error":false},{"allele":0,"quality":{"quality":17},"error":false},{"allele":0,"quality":{"quality":18},"error":false},{"allele":0,"quality":{"quality":25},"error":false}
```
The last field may be less useful for plotting, but describes the observations in JSON format so the specific example can be re-run if desired.

## License

MIT
