var Genotype = require('genotype').Genotype
var Observation = require('genotype').Observation
var multichoose = require('multichoose')

function makeObservations(countQuals) {
  return countQuals.map(function(x) { return new Observation(x[0], x[1]) })
}

var observations =
  makeObservations(
    [ [1, 30],
      [0, 30],
      [1, 10],
      [1, 20] ])

function allGenotypesOfPloidy(ploidy, alleles) {
  var genotypes = []
  multichoose(ploidy, alleles, function(x) {
    genotypes.push(new Genotype(x)) })
  return genotypes
}

var genotypes = allGenotypesOfPloidy(2, [0, 1])

//console.log(genotypes)
//console.log(observations)

genotypes.forEach(function(genotype) {
  console.log(genotype.toString(), genotype.likelihood(observations)) })
genotypes.forEach(function(genotype) {
  console.log(genotype.toString(), genotype.samtoolsLikelihood(observations)) })
