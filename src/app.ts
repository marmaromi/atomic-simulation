import {rVolumeModel, spheresFinalResultsModel} from './models';
import {Spheres} from './spheres';
import XLSX from 'xlsx';
import fs from 'fs';

function runSpheresSimulation(nSpheres: number, rVolumeObj: rVolumeModel, nCollisions: number) {
  let results: spheresFinalResultsModel[] = [];
  for (let rVolume = rVolumeObj.start; rVolume <= rVolumeObj.end; rVolume += rVolumeObj.step) {
    const spheres = new Spheres(nSpheres, rVolume, nCollisions);
    results.push({
      'rVol': spheres.rVolume.toFixed(2),
      'sigma': spheres.sigma.toString(),
      'sum(time)': spheres.totalTCOL.toString(),
      'sum(dv)': spheres.totalDV.toString(),
      'pv0/kT': spheres.pv0.toString(),
    });
  }
  return results;
}

function writeFinalResults(results: spheresFinalResultsModel[], nSpheres: number, nCollisions: number) {
  console.table(results);
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.json_to_sheet(results);
  XLSX.utils.book_append_sheet(wb, ws, 'Table Data');
  fs.writeFileSync(`./output/final_results_${nSpheres}_${nCollisions}.xlsx`, XLSX.write(wb, {type: 'buffer'}));
  console.log('Table data written to file');
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const t0 = performance.now();
const nSpheres = 256;
const nCollisions = 10000;
const rVolume: rVolumeModel = {start: 1.25, end: 1.4, step: 0.05};

const finalResults = runSpheresSimulation(nSpheres, rVolume, nCollisions);
writeFinalResults(finalResults, nSpheres, nCollisions);

const t1 = performance.now();
console.log('Calculations took ' + ((t1 - t0) / 1000).toFixed(2) + ' seconds.');