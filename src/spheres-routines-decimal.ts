import {Calc} from './Calc';
import {spheresResultsPropertiesModel} from './models';
import {Decimal} from 'decimal.js';

const spheresDecimal = (nSpheresNumber: number, rVolumeNumber: number, nCollisionsNumber: number) => {
  Decimal.set({precision: 16, rounding: 4});
  let nSpheres = new Decimal(nSpheresNumber);
  let rVolume = new Decimal(rVolumeNumber);
  let nCollisions = new Decimal(nCollisionsNumber);
  validateInput(nSpheres, rVolume);
  const sigma: Decimal = computeDiameter(nSpheres, rVolume);
  let positions: Decimal[][] = assignPositions(nSpheres);
  let velocities: Decimal[][] = assignVelocitiesZeroMom(nSpheres);
  // writeInitial(nSpheres, rVolume, nCollisions);
  let collisionTime = initCollisionsTable(nSpheres, positions, velocities, sigma);
  let deltaVel: Decimal[] = [];
  let momentum: Decimal[];
  let kineticEnergy: Decimal;
  let vSample: Decimal;
  let orderP: Decimal;
  let dv: Decimal;
  const results: spheresResultsPropertiesModel[] = [];
  let totalDV = new Decimal(0);
  let totalTCOL = new Decimal(0);

  //start collision loop
  /*
  for (let i = 0; i < nCollisionsNumber; i++) {
    let {tCol, iCol, jCol} = retrieveCollisionsInfo(collisionTime);
    ({positions, velocities, deltaVel} = advanceSimulation(nSpheres, positions, velocities, iCol, jCol, tCol, sigma));
    ({momentum, kineticEnergy, vSample, orderP, dv} = computeProperties(positions, velocities, deltaVel, nSpheres));
    results.push(writeProperties(i + 1, momentum, kineticEnergy, tCol, vSample, orderP, dv));
    collisionTime = updateCollisionsTable(nSpheres, positions, velocities, sigma, iCol, jCol, tCol, collisionTime);

    if (i > 5000) { // skip first n collisions
      totalDV = totalDV.plus(dv);
      totalTCOL = totalTCOL.plus(tCol);
    }
  }
  */
  const pv0 = rVolume.mul(totalDV.mul(sigma.div(totalTCOL.mul(nSpheres.mul(3))).plus(1)));
  // writeResults(results, nSpheres, rVolume, nCollisions);
  return {sigma, totalTCOL, totalDV, pv0};
};

const validateInput = (nSpheres: Decimal, rVolume: Decimal) => {
  const nInt = nSpheres.div(4).pow(Decimal.div(1, 3));
  const nIntRounded = Decimal.round(nInt.times(Decimal.pow(10, 12))).div(Decimal.pow(10, 12));
  if (!nIntRounded.isInteger()) {
    throw Error(`Wrong nSpheres!, nSpheres must be a perfect cube of 4. nSpheres = ${nSpheres}`);
  }
  if (rVolume.lt(1)) {
    throw Error(`Wrong rVolume!, rVolume must be >= 1. rVolume = ${rVolume}`);
  }
};

const computeDiameter = (nSpheres: Decimal, rVolume: Decimal, boxLength: Decimal = new Decimal(1)): Decimal => {
  return Decimal.sqrt(2).div(nSpheres.mul(rVolume)).pow(Decimal.div(1, 3)).mul(boxLength);
};

const assignPositions = (nSpheres: Decimal): Decimal[][] => {
  const nInt = nSpheres.div(4).pow(Decimal.div(1, 3)).toNumber();
  const a = new Decimal(1).div(nInt);
  const halfA = a.div(2);

  const positions: Decimal[][] = [];
  for (let i = 0; i < nInt; i++) {
    const ai = a.times(i);
    for (let j = 0; j < nInt; j++) {
      const aj = a.times(j);
      for (let k = 0; k < nInt; k++) {
        const ak = a.times(k);
        positions.push([ai, aj, ak]);
        positions.push([ai.plus(halfA), aj.plus(halfA), ak]);
        positions.push([ai.plus(halfA), aj, ak.plus(halfA)]);
        positions.push([ai, aj.plus(halfA), ak.plus(halfA)]);
      }
    }
  }
  return positions;
};

const vu = [
  {u: 0.67663645744323730, v: 0.23101162910461426},
  {u: 0.61373567581176758, v: 0.055805444717407227},
  {u: 0.92427754402160645, v: 0.33541297912597656},
  {u: 0.28933930397033691, v: 0.92796134948730469},
];

const assignVelocitiesZeroMom = (nSpheres: Decimal): Decimal[][] => {
  const speed = Decimal.sqrt(3);
  const twoPi = Decimal.mul(2, Decimal.acos(-1));
  const velocities: Decimal[][] = [];
  let vSumX = new Decimal(0);
  let vSumY = new Decimal(0);
  let vSumZ = new Decimal(0);
  for (let i = new Decimal(0); i.lt(nSpheres); i = i.plus(1)) {
    const u = new Decimal(vu[i.toNumber()].u) //Decimal.random();
    const v = new Decimal(vu[i.toNumber()].v) //Decimal.random();
    const theta = twoPi.mul(u);
    const phi = Decimal.acos(Decimal.mul(2, v).minus(1));
    const velX = speed.mul(Decimal.sin(phi)).mul(Decimal.cos(theta));
    const velY = speed.mul(Decimal.sin(phi)).mul(Decimal.sin(theta));
    const velZ = speed.mul(Decimal.cos(phi));
    velocities.push([velX, velY, velZ]);
    vSumX = vSumX.plus(velX);
    vSumY = vSumY.plus(velY);
    vSumZ = vSumZ.plus(velZ);
  }

  const vAvgX = vSumX.div(nSpheres);
  const vAvgY = vSumY.div(nSpheres);
  const vAvgZ = vSumZ.div(nSpheres);

  for (let i = 0; i < nSpheres.toNumber(); i++) {
    velocities[i][0] = velocities[i][0].minus(vAvgX);
    velocities[i][1] = velocities[i][1].minus(vAvgY);
    velocities[i][2] = velocities[i][2].minus(vAvgZ);
  }
  return velocities;
};

const writeInitial = (nSpheres: Decimal, rVolume: Decimal, nCollisions: Decimal) => {
  console.group('Input received:');
  console.log('Number of spheres: ', nSpheres.toString());
  console.log('Reduced volume: ', rVolume.toFixed(2));
  console.log('Number of collisions: ', nCollisions.toString());
  // console.log('Initial positions of the spheres: ', this.positions);
  // console.log('Initial velocities of the spheres: ', this.velocities);
  console.groupEnd();
  console.log('----------------------------------------');
};

const initCollisionsTable = (nSpheres: Decimal, positions: Decimal[][], velocities: Decimal[][], sigma: Decimal): Decimal[][] => {
  let collisionTime: Decimal[][] = new Array(nSpheres.toNumber());
  for (let index = 0; index < nSpheres.toNumber(); index++) {
    collisionTime[index] = new Array(nSpheres.toNumber()).fill(new Decimal(Infinity));
  }

  const nSpheresNumber = nSpheres.toNumber();
  for (let i = 1; i <= nSpheresNumber; i++) {
    for (let j = i + 1; j <= nSpheresNumber; j++) {
      const k = i - 1;
      const l = j - 1;
      collisionTime = loopCollisionsTable(k, l, positions, velocities, sigma, collisionTime);
    }
  }
  return collisionTime;
};

const retrieveCollisionsInfo = (collisionTime: Decimal[][]): { tCol: Decimal, iCol: number, jCol: number } => {
  const n = collisionTime.length;
  let tCol = new Decimal(Infinity);
  let iCol = -1;
  let jCol = -1;

  for (let i = 0; i < n - 1; i++) {
    for (let j = i + 1; j < n; j++) {
      if (collisionTime[i][j].lt(tCol)) {
        tCol = collisionTime[i][j];
        iCol = i;
        jCol = j;
      }
    }
  }
  return { tCol, iCol, jCol };
};

const advanceSimulation = (nSpheres: Decimal, positions: Decimal[][], velocities: Decimal[][], iCol: number, jCol: number, tCol: Decimal, sigma: Decimal) => {
  for (let i = new Decimal(0); i.lt(nSpheres); i = i.plus(1)) {
    //update positions
    positions[i.toNumber()] = [
      positions[i.toNumber()][0].plus(velocities[i.toNumber()][0].times(tCol)),
      positions[i.toNumber()][1].plus(velocities[i.toNumber()][1].times(tCol)),
      positions[i.toNumber()][2].plus(velocities[i.toNumber()][2].times(tCol)),
    ];

    // If the new position is outside the simulation wall
    positions[i.toNumber()][0] = positions[i.toNumber()][0].mod(1);
    positions[i.toNumber()][1] = positions[i.toNumber()][1].mod(1);
    positions[i.toNumber()][2] = positions[i.toNumber()][2].mod(1);
  }

  //update velocities of next collision
  let uij = [
    velocities[iCol][0].minus(velocities[jCol][0]),
    velocities[iCol][1].minus(velocities[jCol][1]),
    velocities[iCol][2].minus(velocities[jCol][2]),
  ];
  let minDistance = new Decimal(Infinity);
  let bij: Decimal;
  let rij: Decimal[] = [];

  for (let ix = -1; ix <= 1; ix++) {
    for (let iy = -1; iy <= 1; iy++) {
      for (let iz = -1; iz <= 1; iz++) {
        const translate = [new Decimal(ix), new Decimal(iy), new Decimal(iz)];
        const imaginaryPosition = [
          positions[jCol][0].plus(translate[0]),
          positions[jCol][1].plus(translate[1]),
          positions[jCol][2].plus(translate[2]),
        ];
        const rijImaginary = [
          positions[iCol][0].minus(imaginaryPosition[0]),
          positions[iCol][1].minus(imaginaryPosition[1]),
          positions[iCol][2].minus(imaginaryPosition[2]),
        ];
        const distance = Calc.dotProduct(rijImaginary, rijImaginary);
        if (distance.lt(minDistance)) {
          rij = rijImaginary;
          bij = Calc.dotProduct(rijImaginary, uij);
          minDistance = distance;
        }
      }
    }
  }
  const sigmaSquared = sigma.pow(2);
  const deltaVel: Decimal[] = [
    (bij.div(sigmaSquared)).mul(rij[0]).neg(),
    (bij.div(sigmaSquared)).mul(rij[1]).neg(),
    (bij.div(sigmaSquared)).mul(rij[2]).neg(),
  ];
  velocities[iCol] = [
    velocities[iCol][0].plus(deltaVel[0]),
    velocities[iCol][1].plus(deltaVel[1]),
    velocities[iCol][2].plus(deltaVel[2]),
  ];
  velocities[jCol] = [
    velocities[jCol][0].minus(deltaVel[0]),
    velocities[jCol][1].minus(deltaVel[1]),
    velocities[jCol][2].minus(deltaVel[2]),
  ];

  return {positions, velocities, deltaVel};
};

const computeProperties = (positions: Decimal[][], velocities: Decimal[][], deltaVel: Decimal[], nSpheres: Decimal, boxLength: Decimal = new Decimal(1)) => {
  let momentum: Decimal[] = [new Decimal(0), new Decimal(0), new Decimal(0)];
  let kineticEnergy: Decimal = new Decimal(0);
  let vSample: Decimal = new Decimal(0);
  let orderP: Decimal = new Decimal(0);

  const a = boxLength.div(nSpheres.div(4).pow(Decimal.div(1, 3)));

  for (let i = new Decimal(0); i.lt(nSpheres); i = i.plus(1)) {
    momentum = [
      momentum[0].add(velocities[i.toNumber()][0]),
      momentum[1].add(velocities[i.toNumber()][1]),
      momentum[2].add(velocities[i.toNumber()][2]),
    ];
    const v = Decimal.sqrt(Calc.dotProduct(velocities[i.toNumber()], velocities[i.toNumber()]));
    kineticEnergy = kineticEnergy.add(Decimal.mul(0.5, v.pow(2)));

    if (v.greaterThan(1) && v.lessThan(2)) {
      vSample = vSample.add(1);
    }

    const pi = Decimal.acos(-1);

    orderP = orderP
      .add(Decimal.cos(Decimal.mul(4, pi).mul(positions[i.toNumber()][0]).div(a)))
      .add(Decimal.cos(Decimal.mul(4, pi).mul(positions[i.toNumber()][1]).div(a)))
      .add(Decimal.cos(Decimal.mul(4, pi).mul(positions[i.toNumber()][2]).div(a)));
  }

  momentum = [momentum[0].div(nSpheres), momentum[1].div(nSpheres), momentum[2].div(nSpheres)];
  kineticEnergy = kineticEnergy.div(nSpheres);
  orderP = orderP.div(nSpheres.mul(3));
  const dv = Decimal.sqrt(Calc.dotProduct(deltaVel, deltaVel));

  return {momentum, kineticEnergy, vSample, orderP, dv};
};

const writeProperties = (nthCollision: number, momentum: Decimal[], kineticEnergy: Decimal, tCol: Decimal, vSample: Decimal, orderP: Decimal, dv: Decimal) => {
  return {
    n: nthCollision,
    momX: momentum[0].toFixed(3),
    momY: momentum[1].toFixed(3),
    momZ: momentum[2].toFixed(3),
    kinE: kineticEnergy.toFixed(3),
    tCol: tCol.toNumber(),
    vSample: vSample.toNumber(),
    orderP: orderP.toFixed(3),
    dv: dv.toNumber(),
  };
};

const updateCollisionsTable = (nSpheres: Decimal, positions: Decimal[][], velocities: Decimal[][], sigma: Decimal, iCol: number, jCol: number, tCol: Decimal, collisionTime: Decimal[][]): Decimal[][] => {
  for (let i = 1; i <= nSpheres.toNumber() - 1; i++) {
    for (let j = i + 1; j <= nSpheres.toNumber(); j++) {
      const k = i - 1;
      const l = j - 1;

      collisionTime[k][l] = collisionTime[k][l].minus(tCol);
      if (k === iCol || k === jCol || l === iCol || l === jCol) {
        collisionTime[k][l] = new Decimal(Infinity);
        collisionTime = loopCollisionsTable(k, l, positions, velocities, sigma, collisionTime);
      }
    }
  }
  return collisionTime;
};

const loopCollisionsTable = (i: number, j: number, positions: Decimal[][], velocities: Decimal[][], sigma: Decimal, collisionTime: Decimal[][]): Decimal[][] => {
  const velocitiesI = velocities[i];
  const velocitiesJ = velocities[j];
  const positionsI = positions[i];
  const positionsJ = positions[j];
  const sigmaSquared = sigma.pow(2);
  const collisionTimeIJ = collisionTime[i][j];

  const uijX = velocitiesI[0].minus(velocitiesJ[0]);
  const uijY = velocitiesI[1].minus(velocitiesJ[1]);
  const uijZ = velocitiesI[2].minus(velocitiesJ[2]);
  const uij = [uijX, uijY, uijZ];
  // console.log('uij', uij.toString());
  const uij2 = Calc.dotProduct3(uij, uij);
  // console.log('uij2', uij2.toString());
  for (let ix = -1; ix <= 1; ix++) {
    const translateX = new Decimal(ix).plus(positionsJ[0]);
    const rijX = positionsI[0].minus(translateX);

    for (let iy = -1; iy <= 1; iy++) {
      const translateY = new Decimal(iy).plus(positionsJ[1]);
      const rijY = positionsI[1].minus(translateY);

      for (let iz = -1; iz <= 1; iz++) {
        const translateZ = new Decimal(iz).plus(positionsJ[2]);
        const rijZ = positionsI[2].minus(translateZ);
        const rij = [rijX, rijY, rijZ];
        // console.log('rij', rij.toString());
        const bij = Calc.dotProduct3(rij, uij);
        const cij = Calc.dotProduct3(rij, rij).minus(sigmaSquared);
        // console.log('bij', bij.toString());
        // console.log('cij', cij.toString());
        if (bij.lt(0)) {
          const discriminant = bij.pow(2).minus(uij2.mul(cij));
          // console.log('discriminant', discriminant.toString());

          if (discriminant.gt(0)) {
            const sqrtDiscriminant = discriminant.sqrt();
            const time = new Decimal(-bij).minus(sqrtDiscriminant).div(uij2);
            //console.log('time', time.toString());

            if (time.lt(collisionTimeIJ)) {
              collisionTime[i][j] = time;
            }
          }
        }
      }
    }
  }
  console.log('collisionTime', collisionTime.toString());
  return collisionTime;
};

const writeResults = (results: spheresResultsPropertiesModel[], nSpheres: Decimal, rVolume: Decimal, nCollisions: Decimal) => {
  const fs = require('fs');
  const XLSX = require('xlsx');
  const wb = XLSX.utils.book_new();
  const ws = XLSX.utils.json_to_sheet(results);
  XLSX.utils.book_append_sheet(wb, ws, 'Table Data');
  const date = new Date();
  fs.writeFileSync(`./output/results_${nSpheres.toNumber()}_${rVolume.toFixed(2)}_${nCollisions.toNumber()}_${date.toDateString()}.xlsx`, XLSX.write(wb, {type: 'buffer'}));
  console.log('Table data written to file');
};

export {spheresDecimal};
