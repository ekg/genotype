var multinomialPmf = require('multinomial-pmf')
var multinomialPmfLn = require('multinomial-pmf').log
var multichoose = require('multichoose')
var Phred = require('phred')
Array.prototype.sample = require('array-sample')

function Genotype(a) {
  this.alleles = a.sort()
  this.ploidy = this.alleles.length
  var count = {}
  this.alleles.forEach(function(i) {
    if (isNaN(count[i])) count[i] = 1
    else count[i] += 1
  })
  this.alleleCount = count
  if (Object.keys(count).length === 1) {
    this.homozygous = true
  } else {
    this.homozygous = false
  }
}

Genotype.prototype.hasAllele = function(a) {
  return this.alleleCount[a] > 0
}

Genotype.prototype.toString = function() {
  return this.alleles.join('/')
}

Genotype.prototype.alleleProbs = function() {
  var that = this
  return Object.keys(this.alleleCount).map(function(a) {
    return that.alleleCount[a] / that.ploidy
  })
}

Genotype.prototype.orderedObservationCounts = function(observations) {
  var obsCount = byAlleleObservationCount(this.alleles, observations)
  var counts = []
  return Object.keys(this.alleleCount).map(function(a) {
    return obsCount[a]
  })
}

// an allele and a quality, which is assumed to be log-space p(err) (e.g. phred format)
function Observation(allele, quality) {
  this.allele = allele
  this.quality = new Phred().fromQuality(quality)
  this.error = undefined
}

Observation.prototype.induceQualityDependentError = function(alternateAllele) {
  if (Math.random() < this.quality.toProb()) {
    //this.allele = alternateAllele
    // make a random allele, not just 0/1
    this.allele = [0,1,2,3].sample(1)
    this.error = true
  } else {
    this.error = false
  }
  return this
}

function byAlleleObservationCount(alleles, observations) {
  var counts = {}
  alleles.forEach(function(a) { counts[a] = 0 })
  observations.forEach(function(obs) {
    counts[obs.allele] += 1
  })
  return counts
}

function byAlleleObservationQsum(alleles, observations) {
  var qsums = {}
  alleles.forEach(function(a) { qsums[a] = 0 })
  observations.forEach(function(obs) {
    qsums[obs.allele] += obs.quality.toLog10Prob()
  })
  return qsums
}

function log10(val) {
  return Math.log(val) / Math.LN10
}

function ln2log10(val) {
  return val * Math.LOG10E
}

Genotype.prototype.samplingProbLog10 = function(observations) {
  return ln2log10(
    multinomialPmfLn(
      this.alleleProbs(),
      this.orderedObservationCounts(observations)))
}

Genotype.prototype.orderedSamplingProbLog10 = function(observations) {
  var probs = this.alleleProbs()
  var obsc = this.orderedObservationCounts(observations)
  var lnprob = 0
  for (var i = 0; i < probs.length; ++i) {
    if (obsc[i] > 0) {
      lnprob += log10(probs[i]) * obsc[i]
    }
  }
  return lnprob
}

Genotype.prototype.likelihood = function(alleles, observations) {
  var count = byAlleleObservationCount(alleles, observations)
  var qsums = byAlleleObservationQsum(alleles, observations)
  var qsumOut = 0
  var that = this
  Object.keys(qsums).forEach(function(a) {
    if (!that.hasAllele(a)) {
      qsumOut += qsums[a]
    }
  })
  if (this.homozygous) {
    return qsumOut
  } else {
    return this.samplingProbLog10(observations) + qsumOut
  }
}

Genotype.prototype.samtoolsLikelihood = function(alleles, observations) {
  
  var count = byAlleleObservationCount(alleles, observations)
  var qsums = byAlleleObservationQsum(alleles, observations)
  var qsumOut = 0
  var that = this
  Object.keys(qsums).forEach(function(a) {
    if (!that.hasAllele(a)) {
      qsumOut += qsums[a]
    }
  })
  if (this.homozygous) {
    return qsumOut
  } else {
    return this.orderedSamplingProbLog10(observations) + qsumOut
  }
}

Genotype.prototype.hetDownweightLikelihood = function(alleles, observations, weight) {
  
  var count = byAlleleObservationCount(alleles, observations)
  var qsums = byAlleleObservationQsum(alleles, observations)
  var qsumOut = 0
  var that = this
  Object.keys(qsums).forEach(function(a) {
    if (!that.hasAllele(a)) {
      qsumOut += qsums[a]
    }
  })
  if (this.homozygous) {
    return qsumOut
  } else {
    return this.orderedSamplingProbLog10(observations) * weight + qsumOut
  }
}

module.exports.Genotype = Genotype
module.exports.Observation = Observation
