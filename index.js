var multinomialPmf = require('multinomial-pmf')
var multichoose = require('multichoose')

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
  var obsCount = byAlleleObservationCount(observations)
  var counts = []
  return Object.keys(this.alleleCount).map(function(a) {
    return obsCount[a]
  })
}

// an allele and a quality, which is assumed to be log-space p(err) (e.g. phred format)
function Observation(allele, quality) {
  this.allele = allele
  this.quality = phred2log10(quality)
}

function byAlleleObservationCount(observations) {
  var count = {}
  observations.forEach(function(obs) {
    if (isNaN(count[obs.allele])) count[obs.allele] = 1
    else count[obs.allele] += 1
  })
  return count
}

function byAlleleObservationQsum(observations) {
  var qsums = {}
  observations.forEach(function(obs) {
    if (isNaN(qsums[obs.allele])) qsums[obs.allele] = obs.quality
    else qsums[obs.allele] += 1
  })
  return qsums
}

function log10(val) {
  return Math.log(val) / Math.LN10;
}

function phred2log10(qual) {
  return Math.LN10 * qual * -.1;
}

Genotype.prototype.samplingProbLog10 = function(observations) {
  return log10(
    multinomialPmf(
      this.alleleProbs(),
      this.orderedObservationCounts(observations)))
}

Genotype.prototype.orderedSamplingProbLog10 = function(observations) {
  var probs = this.alleleProbs()
  var obsc = this.orderedObservationCounts(observations)
  var lnprob = 0
  for (var i = 0; i < probs.length; ++i) {
    lnprob += log10(probs[i]) * log10(obsc[i])
  }
  return lnprob
}

Genotype.prototype.likelihood = function(observations) {
  
  var count = byAlleleObservationCount(observations)
  var qsums = byAlleleObservationQsum(observations)
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

Genotype.prototype.samtoolsLikelihood = function(observations) {
  
  var count = byAlleleObservationCount(observations)
  var qsums = byAlleleObservationQsum(observations)
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

module.exports.Genotype = Genotype
module.exports.Observation = Observation
