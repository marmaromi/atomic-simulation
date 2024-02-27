import {Calc} from './Calc';
import {spheresResultsPropertiesModel} from './models';

const spheres = (nSpheres: number, rVolume: number, nCollisions: number) => {
  validateInput(nSpheres, rVolume);
  const sigma: number = computeDiameter(nSpheres, rVolume);
  let positions: number[][] = assignPositions(nSpheres);
  let velocities: number[][] = assignVelocitiesZeroMom(nSpheres);
  writeInitial(nSpheres, rVolume, nCollisions);
  let collisionTime = initCollisionsTable(nSpheres, positions, velocities, sigma);
  let deltaVel: number[] = [];
  let momentum: number[];
  let kineticEnergy: number;
  let vSample: number;
  let orderP: number;
  let dv: number;
  const results: spheresResultsPropertiesModel[] = [];
  let totalDV = 0;
  let totalTCOL = 0;

  //start collision loop
  for (let i = 0; i < nCollisions; i++) {
    let {tCol, iCol, jCol} = retrieveCollisionsInfo(collisionTime);
    ({positions, velocities, deltaVel} = advanceSimulation(nSpheres, positions, velocities, iCol, jCol, tCol, sigma));
    ({momentum, kineticEnergy, vSample, orderP, dv} = computeProperties(positions, velocities, deltaVel, nSpheres));
    results.push(writeProperties(i + 1, momentum, kineticEnergy, tCol, vSample, orderP, dv));
    collisionTime = updateCollisionsTable(nSpheres, positions, velocities, sigma, iCol, jCol, tCol, collisionTime);

    if (i > 5000) { // skip first n collisions
      totalDV += dv;
      totalTCOL += tCol;
    }
  }
  const pv0 = rVolume * (totalDV * (sigma / (totalTCOL * nSpheres * 3)) + 1);
  writeResults(results, nSpheres, rVolume, nCollisions);
  return {sigma, totalTCOL, totalDV, pv0};
};

const validateInput = (nSpheres: number, rVolume: number) => {
  const nInt = Math.pow(nSpheres / 4, 1 / 3);
  const nIntRounded = Math.round(nInt * 1e12) / 1e12;
  if (!Number.isInteger(nIntRounded)) {
    throw Error(`Wrong nSpheres!, nSpheres must be a perfect cube of 4. nSpheres = ${nSpheres}`);
  }
  if (rVolume < 1) {
    throw Error(`Wrong rVolume!, rVolume must be >= 1. rVolume = ${rVolume}`);
  }
};

const computeDiameter = (nSpheres: number, rVolume: number, boxLength: number = 1): number => {
  return Math.pow(Math.sqrt(2) / (nSpheres * rVolume), 1 / 3) * boxLength;
};

const assignPositions = (nSpheres: number): number[][] => {
  const nInt = Math.pow(nSpheres / 4, 1 / 3);
  const a = 1 / nInt;
  const halfA = a / 2;

  const positions = [];
  for (let i = 0; i < nInt; i++) {
    const ai = a * i;
    for (let j = 0; j < nInt; j++) {
      const aj = a * j;
      for (let k = 0; k < nInt; k++) {
        const ak = a * k;
        positions.push([ai, aj, ak]);
        positions.push([ai + halfA, aj + halfA, ak]);
        positions.push([ai + halfA, aj, ak + halfA]);
        positions.push([ai, aj + halfA, ak + halfA]);
      }
    }
  }
  return positions;
};

const assignVelocitiesZeroMom = (nSpheres: number): number[][] => {
  const speed = Math.sqrt(3);
  const twoPi = 2 * Math.PI;
  const velocities = [];
  let vSumX = 0;
  let vSumY = 0;
  let vSumZ = 0;
  for (let i = 0; i < nSpheres; i++) {
    const u = Math.random();
    const v = Math.random();
    const theta = twoPi * u;
    const phi = Math.acos(2 * v - 1);
    const velX = speed * Math.sin(phi) * Math.cos(theta);
    const velY = speed * Math.sin(phi) * Math.sin(theta);
    const velZ = speed * Math.cos(phi);
    velocities.push([velX, velY, velZ]);
    vSumX += velX;
    vSumY += velY;
    vSumZ += velZ;
  }
  const vAvgX = vSumX / nSpheres;
  const vAvgY = vSumY / nSpheres;
  const vAvgZ = vSumZ / nSpheres;

  for (let i = 0; i < nSpheres; i++) {
    velocities[i][0] -= vAvgX;
    velocities[i][1] -= vAvgY;
    velocities[i][2] -= vAvgZ;
  }
  return velocities;
};

const writeInitial = (nSpheres: number, rVolume: number, nCollisions: number) => {
  console.group('Input received:');
  console.log('Number of spheres: ', nSpheres);
  console.log('Reduced volume: ', rVolume.toFixed(2));
  console.log('Number of collisions: ', nCollisions);
  // console.log('Initial positions of the spheres: ', this.positions);
  // console.log('Initial velocities of the spheres: ', this.velocities);
  console.groupEnd();
  console.log('----------------------------------------');
};

const initCollisionsTable = (nSpheres: number, positions: number[][], velocities: number[][], sigma: number): number[][] => {
  let collisionTime: number[][] = new Array(nSpheres);
  for (let index = 0; index < nSpheres; index++) {
    collisionTime[index] = new Array(nSpheres).fill(Infinity);
  }

  for (let i = 1; i <= nSpheres - 1; i++) {
    for (let j = i + 1; j <= nSpheres; j++) {
      const k = i - 1;
      const l = j - 1;
      collisionTime = loopCollisionsTable(k, l, positions, velocities, sigma, collisionTime);
    }
  }
  return collisionTime;
};

const retrieveCollisionsInfo = (collisionTime: number[][]) => {
  const n = collisionTime.length;
  let tCol = Infinity;
  let iCol = -1;
  let jCol = -1;

  for (let i = 0; i < n - 1; i++) {
    for (let j = i + 1; j < n; j++) {
      if (collisionTime[i][j] < tCol) {
        tCol = collisionTime[i][j];
        iCol = i;
        jCol = j;
      }
    }
  }
  return {tCol, iCol, jCol};
};

const advanceSimulation = (nSpheres: number, positions: number[][], velocities: number[][], iCol: number, jCol: number, tCol: number, sigma: number) => {
  for (let i = 0; i < nSpheres; i++) {
    //update positions
    positions[i] = [
      positions[i][0] + velocities[i][0] * tCol,
      positions[i][1] + velocities[i][1] * tCol,
      positions[i][2] + velocities[i][2] * tCol,
    ];

    // If the new position is outside the simulation wall
    positions[i][0] = positions[i][0] % 1;
    positions[i][1] = positions[i][1] % 1;
    positions[i][2] = positions[i][2] % 1;
  }

  //update velocities of next collision
  const uij = [
    velocities[iCol][0] - velocities[jCol][0],
    velocities[iCol][1] - velocities[jCol][1],
    velocities[iCol][2] - velocities[jCol][2],
  ];
  let minDistance = Infinity;
  let bij;
  let rij: number[] = [];

  for (let ix = -1; ix <= 1; ix++) {
    for (let iy = -1; iy <= 1; iy++) {
      for (let iz = -1; iz <= 1; iz++) {
        const translate = [ix, iy, iz];
        const imaginaryPosition = [
          positions[jCol][0] + translate[0],
          positions[jCol][1] + translate[1],
          positions[jCol][2] + translate[2],
        ];
        const rijImaginary = [
          positions[iCol][0] - imaginaryPosition[0],
          positions[iCol][1] - imaginaryPosition[1],
          positions[iCol][2] - imaginaryPosition[2],
        ];
        const distance = Calc.dotProduct2(rijImaginary, rijImaginary);
        if (distance < minDistance) {
          rij = rijImaginary;
          bij = Calc.dotProduct2(rijImaginary, uij);
          minDistance = distance;
        }
      }
    }
  }
  const sigmaSquared = Math.pow(sigma, 2);
  const deltaVel: number[] = [
    -(bij / sigmaSquared) * rij[0],
    -(bij / sigmaSquared) * rij[1],
    -(bij / sigmaSquared) * rij[2],
  ];
  velocities[iCol] = [
    velocities[iCol][0] + deltaVel[0],
    velocities[iCol][1] + deltaVel[1],
    velocities[iCol][2] + deltaVel[2],
  ];
  velocities[jCol] = [
    velocities[jCol][0] - deltaVel[0],
    velocities[jCol][1] - deltaVel[1],
    velocities[jCol][2] - deltaVel[2],
  ];

  return {positions, velocities, deltaVel};
};

const computeProperties = (positions: number[][], velocities: number[][], deltaVel: number[], nSpheres: number, boxLength = 1) => {
  let momentum = [0, 0, 0];
  let kineticEnergy = 0;
  let vSample = 0;
  let orderP = 0;
  let dv: number;

  const a = boxLength / Math.pow(nSpheres / 4, 1 / 3);

  for (let i = 0; i < nSpheres; i++) {
    momentum = [
      momentum[0] + velocities[i][0],
      momentum[1] + velocities[i][1],
      momentum[2] + velocities[i][2],
    ];
    const v = Math.sqrt(Calc.dotProduct2(velocities[i], velocities[i]));
    kineticEnergy += 0.5 * v ** 2;

    if (1 <= v && v < 2) {
      vSample += 1;
    }

    orderP += Math.cos(4 * Math.PI * positions[i][0] / a) +
      Math.cos(4 * Math.PI * positions[i][1] / a) +
      Math.cos(4 * Math.PI * positions[i][2] / a);
  }

  momentum = [momentum[0] / nSpheres, momentum[1] / nSpheres, momentum[2] / nSpheres];
  kineticEnergy /= nSpheres;
  orderP /= 3 * nSpheres;
  dv = Math.sqrt(Calc.dotProduct2(deltaVel, deltaVel));

  return {momentum, kineticEnergy, vSample, orderP, dv};
};

const writeProperties = (nthCollision: number, momentum: number[], kineticEnergy: number, tCol: number, vSample: number, orderP: number, dv: number) => {
  return {
    n: nthCollision,
    momX: momentum[0].toFixed(3),
    momY: momentum[1].toFixed(3),
    momZ: momentum[2].toFixed(3),
    kinE: kineticEnergy.toFixed(3),
    tCol: tCol,
    vSample: vSample,
    orderP: orderP.toFixed(3),
    dv: dv,
  };
};

const updateCollisionsTable = (nSpheres: number, positions: number[][], velocities: number[][], sigma: number, iCol: number, jCol: number, tCol: number, collisionTime: number[][]) => {
  for (let i = 1; i <= nSpheres - 1; i++) {
    for (let j = i + 1; j <= nSpheres; j++) {
      const k = i - 1;
      const l = j - 1;

      collisionTime[k][l] -= tCol;
      if (k === iCol || k === jCol || l === iCol || l === jCol) {
        collisionTime[k][l] = Infinity;
        collisionTime = loopCollisionsTable(k, l, positions, velocities, sigma, collisionTime);
      }
    }
  }
  return collisionTime;
};

const loopCollisionsTable = (i: number, j: number, positions: number[][], velocities: number[][], sigma: number, collisionTime: number[][]): number[][] => {
  const velocitiesI = velocities[i];
  const velocitiesJ = velocities[j];
  const positionsI = positions[i];
  const positionsJ = positions[j];
  const sigmaSquared = sigma * sigma;
  const collisionTimeIJ = collisionTime[i][j];

  for (let ix = -1; ix <= 1; ix++) {
    const translateX = ix + positionsJ[0];

    for (let iy = -1; iy <= 1; iy++) {
      const translateY = iy + positionsJ[1];

      for (let iz = -1; iz <= 1; iz++) {
        const translateZ = iz + positionsJ[2];
        const rijX = positionsI[0] - translateX;
        const rijY = positionsI[1] - translateY;
        const rijZ = positionsI[2] - translateZ;
        const uijX = velocitiesI[0] - velocitiesJ[0];
        const uijY = velocitiesI[1] - velocitiesJ[1];
        const uijZ = velocitiesI[2] - velocitiesJ[2];
        const bij = rijX * uijX + rijY * uijY + rijZ * uijZ;
        const cij = rijX * rijX + rijY * rijY + rijZ * rijZ - sigmaSquared;

        if (bij < 0) {
          const uij2 = uijX * uijX + uijY * uijY + uijZ * uijZ;
          const discriminant = bij * bij - uij2 * cij;

          if (discriminant > 0) {
            const sqrtDiscriminant = Math.sqrt(discriminant);
            const time = (-bij - sqrtDiscriminant) / uij2;

            if (time < collisionTimeIJ) {
              collisionTime[i][j] = time;
            }
          }
        }
      }
    }
  }
  return collisionTime;
};

const writeResults = (results: spheresResultsPropertiesModel[], nSpheres: number, rVolume: number, nCollisions: number) => {
  const fs = require('fs');
  const XLSX = require('xlsx');
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.json_to_sheet(results);
  XLSX.utils.book_append_sheet(wb, ws, 'Table Data');
  const date = new Date();
  fs.writeFileSync(`./output/results_${nSpheres}_${rVolume.toFixed(2)}_${nCollisions}_${date.toDateString()}.xlsx`, XLSX.write(wb, {type: 'buffer'}));
  console.log('Table data written to file');
};

export {spheres};
