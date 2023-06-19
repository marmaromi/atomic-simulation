import {Spheres} from './spheres';

const t0 = performance.now();
const results = [];
for (let rVolume = 1; rVolume <= 1; rVolume += 0.05) {
  const spheres = new Spheres(4, rVolume, 100);
  results.push({
    'rVol': spheres.rVolume.toFixed(2),
    'sigma': spheres.sigma,
    'sum(time)': spheres.totalTCOL,
    'sum(dv)': spheres.totalDV,
    'pv0/kT': spheres.pv0,
  });
}
console.table(results);
const t1 = performance.now();
console.log('calculations took ' + (t1 - t0) + ' milliseconds.');
