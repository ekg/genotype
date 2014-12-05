# genotype

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

## License

MIT
