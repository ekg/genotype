#!/usr/bin/env node

var Genotype = require('genotype').Genotype
var Observation = require('genotype').Observation
var multichoose = require('multichoose')
var Phred = require('phred')
var range = require('range')

Array.prototype.sample = require('array-sample')

function makeObservations(countQuals) {
  return countQuals.map(function(x) { return new Observation(x[0], x[1]) })
}

function allGenotypesOfPloidy(ploidy, alleles) {
  var genotypes = []
  multichoose(ploidy, alleles, function(x) {
    genotypes.push(new Genotype(x)) })
  return genotypes
}

function simulateObservations(genotype, qualities, number) {
  return genotype.alleles.sample(number).map(function(allele) {
    var qual = qualities.sample(1)
    
  })
}

function normalizeLikelihoods(likelihoods) {
  var maxq = likelihoods.reduce(function(x,y) {
    if (x > y) return x
    else return y
  })
  return likelihoods.map(function(x) { return x - maxq })
}

/*
var qobs =
    [ [1, 30],
      [0, 30],
      [1, 10],
      [1, 20] ]
*/

// where is zip.... sheesh
function genotypesAndLikelihoods(genotypes, likelihoods) {
  var result = []
  for (var i = 0; i < genotypes.length; ++i) {
    result.push(genotypes[i], likelihoods[i])
  }
  return result
}

function freebayesLikelihoods(genotypes, alleles, observations) {
  var likelihoods = normalizeLikelihoods(genotypes.map(function(genotype) { return genotype.likelihood(alleles, observations) }))
  return likelihoods //.map(function(x) { return new Phred().fromLog10Prob(x).toProb() })
}

function samtoolsLikelihoods(genotypes, alleles, observations) {
  var likelihoods = normalizeLikelihoods(genotypes.map(function(genotype) { return genotype.samtoolsLikelihood(alleles, observations) }))
  return likelihoods //.map(function(x) { return new Phred().fromLog10Prob(x).toProb() })
}

if (process.argv.length < 4) {
  console.log('usage:', process.argv[1], '[reps]', '[depth factor]')
  console.log('tests samtools and freebayes genotype likelihoods by simulating observations')
  console.log('outputs a descriptive summary of results')
  process.exit(1)
}

var reps = process.argv[2]
var coverageFactor = process.argv[3]


// write header
console.log(['method',
             'true.genotype',
             'status',
             'called.genotype',
             'true.genotype.likelihood',
             'n.observations',
             'p.homref',
             'p.het',
             'p.homalt',
             'observations'].join('\t'))

function formatTestResult(method, genotypes, trueGenotype, trueGenotypeIndex, observations, likelihoods) {
  var calledGenotypeIndex
  for (var i = 0; i < likelihoods.length; ++i) {
    if (likelihoods[i] === 0) {
      calledGenotypeIndex = i
      break
    }
  }
  return [method,
          trueGenotype.toString(),
          likelihoods[trueGenotypeIndex] === 0 ? "pass" : "fail",
          genotypes[calledGenotypeIndex].toString(),
          likelihoods[trueGenotypeIndex],
          observations.length,
          likelihoods.join('\t'),
          JSON.stringify(observations)].join('\t')
}


range(0,reps).forEach(function() {
  var result = {}
  var ploidy = 2
  var alleles = [0, 1]
  var possibleAlleles = [0, 1, 2, 3]
  var qualities = range(1,30)
  for (var i = 0; i < 10; ++i) qualities.push(30)
  var genotype = new Genotype(alleles.sample(2))
  var numObservations = Math.round(Math.random() * coverageFactor)
  while (numObservations < 2) numObservations = Math.round(Math.random() * coverageFactor)
  var obsQualities = qualities.sample(numObservations) 
  var alleleObservations = genotype.alleles.sample(numObservations)
  var observations = []
  for (var i = 0; i < numObservations; ++i) {
    var obs = new Observation(alleleObservations[i], obsQualities[i])
    var altAllele
    do {
      altAllele = alleles.sample(1)
    } while (altAllele === obs.allele)
    obs.induceQualityDependentError(altAllele)
    observations.push(obs)
  }
  /*
  result.observations = observations
  result.likelihoods = likelihoods
  */
  var genotypeIndex = 0
  var genotypes = allGenotypesOfPloidy(ploidy, alleles)
  for (var i = 0; i < genotypes.length; ++i) {
    if (genotypes[i].toString() === genotype.toString()) {
      genotypeIndex = i
      break
    }
  }
  var genotypes = allGenotypesOfPloidy(ploidy, alleles)
  console.log(formatTestResult('freebayes', genotypes, genotype, genotypeIndex, observations, freebayesLikelihoods(genotypes, possibleAlleles, observations)))
  console.log(formatTestResult('samtools', genotypes, genotype, genotypeIndex, observations, samtoolsLikelihoods(genotypes, possibleAlleles, observations)))
  console.log(formatTestResult('hetdown', genotypes, genotype, genotypeIndex, observations, hetDownweightLikelihoods(genotypes, possibleAlleles, observations)))
})

