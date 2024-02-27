import {rVolumeModel, spheresFinalResultsModel} from './models';
import {Spheres} from './spheres';
import XLSX from 'xlsx';
import fs from 'fs';
import {Spheres2} from './spheres2';
import {Spheres3} from './spheres3';
import {SpheresCopilot} from './spheres-copilot';
import {spheres} from './spheres-routines';
import {spheresDecimal} from './spheres-routines-decimal';
import {Decimal} from 'decimal.js';

function runSpheresSimulation(nSpheres: number, rVolumeObj: rVolumeModel, nCollisions: number) {
  let results: spheresFinalResultsModel[] = [];
  let sigma: Decimal;
  let totalTCOL: Decimal;
  let totalDV: Decimal;
  let pv0: Decimal;
  for (let rVolume = rVolumeObj.start; rVolume <= rVolumeObj.end; rVolume += rVolumeObj.step) {
    const t0 = performance.now();
    // const spheres = new Spheres(nSpheres, rVolume, nCollisions);
    ({sigma, totalTCOL, totalDV, pv0} = spheresDecimal(nSpheres, rVolume, nCollisions));
    // results.push({
    //   'rVol': rVolume.toFixed(2),
    //   'sigma': sigma.toString(),
    //   'sum(time)': totalTCOL.toString(),
    //   'sum(dv)': totalDV.toString(),
    //   'pv0/kT': pv0.toString(),
    // });
    // const t1 = performance.now();
    // console.log('Time to run simulation: ' + ((t1 - t0) / 1000).toFixed(2) + ' seconds.');
  }
  return results;
}

function writeFinalResults(results: spheresFinalResultsModel[], nSpheres: number, nCollisions: number) {
  console.table(results);
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.json_to_sheet(results);
  XLSX.utils.book_append_sheet(wb, ws, 'Table Data');
  const date = new Date();
  fs.writeFileSync(`./output/final_results_${nSpheres}_${nCollisions}_${date.toDateString()}.xlsx`, XLSX.write(wb, {type: 'buffer'}));
  console.log('Table data written to file');
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const t0 = performance.now();
const nSpheres = 4;
const nCollisions = 1;
const rVolume: rVolumeModel = {start: 1, end: 1, step: 0.05};

const finalResults = runSpheresSimulation(nSpheres, rVolume, nCollisions);
writeFinalResults(finalResults, nSpheres, nCollisions);

const t1 = performance.now();
console.log('Total Calculation Time: ' + ((t1 - t0) / 1000).toFixed(2) + ' seconds.');