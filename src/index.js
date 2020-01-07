'use strict';

const fs = require('fs');
const randomNormal = require('random-normal');
const predictor = require('nmr-predictor');
const simulation = require('nmr-simulation');

var defaultOptions = {
  frequency: 400,
  from: 1,
  to: 10,
  lineWidth: 1,
  nbPoints: 1024 * 1,
  maxClusterSize: 7,
  output: 'y',
  withNoise: false
};

async function makeSpinusPredictions(inputPath, outpuPath) {
  let molfiles = fs.readdirSync(inputPath);
  for (let file of molfiles) {
    let molfile = fs.readFileSync(inputPath + "/" + file).toString();
    let prediction = await predictor.spinus(molfile);
    fs.writeFileSync(outpuPath + '/' + file.replace('.mol', '.json'), JSON.stringify(prediction));
    console.log(file);
  }
}

function generateConcentrations(elements, nSamples) {
  let result = new Array(nSamples);
  for (let i = 0; i < nSamples; i++) {
    // It works if the concentrations are independent for each pair of pure compounds
    result[i] = elements.map(element => element.dist.nextRandom())
  }

  return result;
}

function generateDistributions(elements, defaultDist = { mean: 1, dev: 0.1, name: 'normal' }) {
  return elements.map(element => {
    let params = Object.assign({}, defaultDist, element.dist);
    params.nextRandom = () => {
      // Negative values are not allowed in this framework!
      let value = Math.round(randomNormal(params) * 10000) / 10000;
      return value < 0 ? 0 : value;
    }

    return Object.assign({}, element, {dist: params});
  });
}

function loadPredictions(inputPath) {
  let result = [];
  let files = fs.readdirSync(inputPath);
  for (let file of files) {
    result.push({name: file, value: JSON.parse(fs.readFileSync(inputPath + "/" + file).toString())});
  }

  return result;
}

function generateMixtures(elements, concentrations, mutator, options1h, outpuPath) {
  concentrations.forEach((row, rowIndex) => {
    console.log(rowIndex);
    let mixture = row.map((weight, index) => {
      if (weight > 0) {
        //console.log(mutator(elements[index].value))
        let spinSystem = simulation.SpinSystem.fromPrediction(mutator(elements[index].value));
        spinSystem.ensureClusterSize(options1h);
        let spectrum = simulation.simulate1D(spinSystem, options1h);
        let norm = weight / spectrum.reduce((sum, val) => sum + val, 0);
        return spectrum.map(value => value * norm)
      } else {
        let spectrum = new Array(options1h.nbPoints);
        for(let i = 0; i < spectrum.length; i++) {
          spectrum[i] = 0;
        }
        return spectrum;
      }
    }).reduce((sum, value) => {
      if(sum.length === 0) {
        sum = new Array(value.length);
        for(let i = 0; i < sum.length; i++) {
          sum[i] = 0;
        }
      }
      for(let i = 0; i < sum.length; i++) {
        sum[i] += value[i]; 
      }
      return sum;
    }, []);
    fs.writeFileSync(outpuPath + rowIndex + ".json" , JSON.stringify(mixture));
  });
}

//STEP 1
//makeSpinusPredictions('data/molecules/', 'data/predictions/');

//STEP 2
let elements = loadPredictions('data/predictions/');
let elements2 = generateDistributions(elements, { mean: 1, dev: 0.5, name: 'normal' });
let concentrations = generateConcentrations(elements2, 30000);

fs.writeFileSync('data/output/elements.json', JSON.stringify(elements2));
fs.writeFileSync('data/output/concentrations.json', JSON.stringify(concentrations));

//STEP 3
let spectra = generateMixtures(elements2, concentrations, x => x, defaultOptions, 'data/output/mixtures/' );
fs.writeFileSync('data/output/spectra.json', JSON.stringify(spectra));



