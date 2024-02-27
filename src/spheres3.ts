import {Calc} from './Calc';
import {Decimal} from 'decimal.js';
import {Big} from 'big.js';
import {spheresResultsPropertiesModel} from './models';

export class Spheres3 {

  nSpheres: Big;
  rVolume: number;
  nCollisions: number;
  sigma: Big;
  positions: Big[][] = [];
  velocities: Big[][] = [];
  collisionTime: Big[][] = [];
  tCol: Big;
  iCol: number;
  jCol: number;
  deltaVel: Big[];
  momentum: Big[];
  kineticEnergy: Big;
  vSample: Big;
  orderP: Big;
  dv: Big;
  pv0: Big;
  totalDV: Big = new Big(0);
  totalTCOL: Big = new Big(0);
  results: spheresResultsPropertiesModel[] = [];

  boxLength = 1;
  // PI = Big.acos(-1);
  infinity = new Big(1e100);

  constructor(nSpheres: number, rVolume: number, nCollisions: number) {
    this.setInput(nSpheres, rVolume, nCollisions);
    this.validateInput();
    this.computeDiameter();
    this.assignPositions();
    this.assignVelocitiesZeroMom();
    this.writeInitial();
    this.updateCollisionsTable('init');

    const tempCalculation = this.totalTCOL.times(this.nSpheres).times(3);
    let skipCollisions = false;

    // start collision loop
    for (let i = 0; i < this.nCollisions; i++) {
      this.retrieveCollisionsInfo();
      this.advanceSimulation();
      this.computeProperties();
      this.writeProperties(i + 1);
      this.updateCollisionsTable('update');

      if (!skipCollisions) {
        this.totalDV = this.totalDV.plus(this.dv);
        this.totalTCOL = this.totalTCOL.plus(this.tCol);
      }

      if (i === 5000) {
        skipCollisions = true;
      }
    }

    const decimalSigma = new Big(this.sigma);
    const pv0Multiplier = this.totalDV.times(decimalSigma).div(tempCalculation).plus(1);
    this.pv0 = new Big(this.rVolume).times(pv0Multiplier);
    this.writeResults();
  }

  setInput(nSpheres: number, rVolume: number, nCollisions: number) {
    this.nSpheres = new Big(nSpheres);
    this.rVolume = rVolume;
    this.nCollisions = nCollisions;
  }

  validateInput() {
    const nInt = Math.pow(this.nSpheres.toNumber() / 4, 1 / 3);
    const nIntRounded = Math.round(nInt * 1e12) / 1e12;
    if (!Number.isInteger(nIntRounded)) {
      throw Error(`Wrong nSpheres!, nSpheres must be a perfect cube of 4. nSpheres = ${this.nSpheres}`);
    }
    if (this.rVolume < 1) {
      throw Error(`Wrong rVolume!, rVolume must be >= 1. rVolume = ${this.rVolume}`);
    }
  }

  computeDiameter() {
    const numSigma = Math.pow(Math.sqrt(2) / (this.nSpheres.toNumber() * this.rVolume), 1/3) * this.boxLength;
    this.sigma = new Big(numSigma);
  }

  assignPositions() {
    const nInt = Math.pow(this.nSpheres.toNumber() / 4, 1 / 3);
    const a = 1 / nInt;
    const halfA = a / 2;

    this.positions = new Array(nInt * nInt * nInt * 4);

    let index = 0;
    for (let i = 0; i < nInt; i++) {
      const xNum = a * i;
      for (let j = 0; j < nInt; j++) {
        const yNum = a * j;
        for (let k = 0; k < nInt; k++) {
          const zNum = a * k;

          const x = new Big(xNum);
          const y = new Big(yNum);
          const z = new Big(zNum);

          this.positions[index++] = [x, y, z];
          this.positions[index++] = [x.add(halfA), y.add(halfA), z];
          this.positions[index++] = [x.add(halfA), y, z.add(halfA)];
          this.positions[index++] = [x, y.add(halfA), z.add(halfA)];
        }
      }
    }
  }

  assignVelocitiesZeroMom() {
    const speed = Math.sqrt(3);
    const twoPi = 2 * Math.PI;
    const nSpheresNum = this.nSpheres.toNumber();

    this.velocities = new Array(nSpheresNum);

    let vSumX = 0;
    let vSumY = 0;
    let vSumZ = 0;

    for (let i = 0; i < nSpheresNum; i++) {
      const u = Math.random();
      const v = Math.random();
      const theta = twoPi * u;
      const phi = Math.acos(2 * v - 1);

      const vX = speed * Math.sin(phi) * Math.cos(theta);
      const vY = speed * Math.sin(phi) * Math.sin(theta);
      const vZ =  speed * Math.cos(phi);

      this.velocities[i] = [new Big(vX), new Big(vY), new Big(vZ)];

      vSumX += vX;
      vSumY += vY;
      vSumZ += vZ;
    }

    const vAvgX = new Big(vSumX / nSpheresNum);
    const vAvgY = new Big(vSumY / nSpheresNum);
    const vAvgZ = new Big(vSumZ / nSpheresNum);

    for (let i = 0; i < nSpheresNum; i++) {
      this.velocities[i][0] = this.velocities[i][0].minus(vAvgX);
      this.velocities[i][1] = this.velocities[i][1].minus(vAvgY);
      this.velocities[i][2] = this.velocities[i][2].minus(vAvgZ);
    }
  }

  writeInitial() {
    console.group('Input received:');
    console.log('Number of spheres: ', this.nSpheres);
    console.log('Reduced volume: ', this.rVolume.toFixed(2));
    console.log('Number of collisions: ', this.nCollisions);
    // console.log('Initial positions of the spheres: ', this.positions);
    // console.log('Initial velocities of the spheres: ', this.velocities);
    console.groupEnd();
    console.log('----------------------------------------');
  }

  retrieveCollisionsInfo() {
    const n = this.collisionTime.length;
    this.tCol = this.infinity;

    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const collisionTime = this.collisionTime[i][j];
        if (collisionTime < this.tCol) {
          this.tCol = collisionTime;
          this.iCol = i;
          this.jCol = j;
        }
      }
    }
  }

  advanceSimulation() {
    const sigmaSquared = this.sigma.pow(2);
    const nSpheres = this.nSpheres.toNumber();
    const iColPositions = this.positions[this.iCol];
    const jColPositions = this.positions[this.jCol];
    const iColVelocities = this.velocities[this.iCol];
    const jColVelocities = this.velocities[this.jCol];
    const uij = [
      iColVelocities[0].minus(jColVelocities[0]),
      iColVelocities[1].minus(jColVelocities[1]),
      iColVelocities[2].minus(jColVelocities[2]),
    ];

    for (let i = 0; i < nSpheres; i++) {
      // Update positions
      const tColTimesVelocity = [
        this.tCol.times(this.velocities[i][0]),
        this.tCol.times(this.velocities[i][1]),
        this.tCol.times(this.velocities[i][2]),
      ];
      this.positions[i] = [
        this.positions[i][0].add(tColTimesVelocity[0]),
        this.positions[i][1].add(tColTimesVelocity[1]),
        this.positions[i][2].add(tColTimesVelocity[2]),
      ];

      // If the new position is outside the simulation wall
      this.positions[i][0] = this.positions[i][0].mod(1);
      this.positions[i][1] = this.positions[i][1].mod(1);
      this.positions[i][2] = this.positions[i][2].mod(1);
    }

    // Find minimum distance and update velocities of next collision
    let minDistance = this.infinity;
    let bij;
    let rij;

    for (let ix = -1; ix <= 1; ix++) {
      for (let iy = -1; iy <= 1; iy++) {
        for (let iz = -1; iz <= 1; iz++) {
          const translate = [ix, iy, iz];
          const imaginaryPosition = [
            jColPositions[0].add(translate[0]),
            jColPositions[1].add(translate[1]),
            jColPositions[2].add(translate[2]),
          ];
          const rijImaginary = [
            iColPositions[0].minus(imaginaryPosition[0]),
            iColPositions[1].minus(imaginaryPosition[1]),
            iColPositions[2].minus(imaginaryPosition[2]),
          ];
          const distance = Calc.dotProduct4(rijImaginary, rijImaginary);
          if (distance.lt(minDistance)) {
            rij = rijImaginary;
            bij = Calc.dotProduct4(rijImaginary, uij);
            minDistance = distance;
          }
        }
      }
    }

    const factor = bij.div(sigmaSquared).neg();
    this.deltaVel = [
      factor.times(rij[0]),
      factor.times(rij[1]),
      factor.times(rij[2]),
    ];

    iColVelocities[0] = iColVelocities[0].add(this.deltaVel[0]);
    iColVelocities[1] = iColVelocities[1].add(this.deltaVel[1]);
    iColVelocities[2] = iColVelocities[2].add(this.deltaVel[2]);
    jColVelocities[0] = jColVelocities[0].minus(this.deltaVel[0]);
    jColVelocities[1] = jColVelocities[1].minus(this.deltaVel[1]);
    jColVelocities[2] = jColVelocities[2].minus(this.deltaVel[2]);

    // Update velocities in the array
    this.velocities[this.iCol] = iColVelocities;
    this.velocities[this.jCol] = jColVelocities;
  }

  computeProperties() {
    this.momentum = [new Big(0), new Big(0), new Big(0)];
    this.kineticEnergy = new Big(0);
    this.vSample = new Big(0);
    this.orderP = new Big(0);

    const nSpheres = this.nSpheres.toNumber();

    const a = this.boxLength / Math.pow(nSpheres / 4, 1 / 3);


    for (let i = 0; i < nSpheres; i++) {
      const [vx, vy, vz] = this.velocities[i];
      const [x, y, z] = this.positions[i];

      this.momentum = [
        this.momentum[0].add(vx),
        this.momentum[1].add(vy),
        this.momentum[2].add(vz),
      ];

      const v = Calc.dotProduct4(this.velocities[i], this.velocities[i]).sqrt();
      this.kineticEnergy = this.kineticEnergy.add(v.mul(v).div(2));

      if (v.gt(1) && v.lt(2)) {
        this.vSample = this.vSample.add(1);
      }

      const orderPNumber =
        Math.cos(x.mul(4 * Math.PI).div(a).toNumber()) +
        Math.cos(y.mul(4 * Math.PI).div(a).toNumber()) +
        Math.cos(z.mul(4 * Math.PI).div(a).toNumber());

      this.orderP = this.orderP.add(new Big(orderPNumber));
    }

    this.momentum = [this.momentum[0].div(this.nSpheres), this.momentum[1].div(this.nSpheres), this.momentum[2].div(this.nSpheres)];
    this.kineticEnergy = this.kineticEnergy.div(this.nSpheres);
    this.orderP = this.orderP.div(this.nSpheres.mul(3));
    const dotProductDeltaVel = Calc.dotProduct4(this.deltaVel, this.deltaVel)
    this.dv = dotProductDeltaVel.sqrt();
  }

  writeProperties(nthCollision: number) {
    this.results.push({
      n: nthCollision,
      momX: this.momentum[0].toFixed(3),
      momY: this.momentum[1].toFixed(3),
      momZ: this.momentum[2].toFixed(3),
      kinE: this.kineticEnergy.toFixed(3),
      tCol: this.tCol.toNumber(),
      vSample: this.vSample.toNumber(),
      orderP: this.orderP.toFixed(3),
      dv: this.dv.toNumber(),
    });
  }

  updateCollisionsTable(role: string) {
    const nSpheresValue = this.nSpheres.toNumber();
    const infinity = this.infinity;

    if (role === 'init') {
      this.collisionTime = new Array(nSpheresValue);

      for (let i = 0; i < nSpheresValue; i++) {
        this.collisionTime[i] = new Array(nSpheresValue);
        for (let j = 0; j < nSpheresValue; j++) {
          this.collisionTime[i][j] = infinity;
        }
      }
    }

    for (let i = 1; i <= nSpheresValue - 1; i++) {
      const k = i - 1;
      for (let j = i + 1; j <= nSpheresValue; j++) {
        const l = j - 1;
        if (role === 'init') {
          this.loopCollisionsTable(k, l);
        } else if (role === 'update') {
          this.collisionTime[k][l] = this.collisionTime[k][l].minus(this.tCol);
          if (k === this.iCol || k === this.jCol || l === this.iCol || l === this.jCol) {
            this.collisionTime[k][l] = infinity;
            this.loopCollisionsTable(k, l);
          }
        }
      }
    }
  }

  loopCollisionsTable(i: number, j: number) {
    const velocitiesI = this.velocities[i];
    const velocitiesJ = this.velocities[j];
    const positionsI = this.positions[i];
    const positionsJ = this.positions[j];
    const uij = [
      velocitiesI[0].minus(velocitiesJ[0]),
      velocitiesI[1].minus(velocitiesJ[1]),
      velocitiesI[2].minus(velocitiesJ[2]),
    ];

    const sigmaSquared = this.sigma.pow(2);

    let uij2 = new Big(0);
    for (let k = 0; k < 3; k++) {
      uij2 = uij2.plus(uij[k].times(uij[k]));
    }

    for (let ix = -1; ix <= 1; ix++) {
      for (let iy = -1; iy <= 1; iy++) {
        for (let iz = -1; iz <= 1; iz++) {
          const imaginaryPositionX = positionsJ[0].plus(ix);
          const imaginaryPositionY = positionsJ[1].plus(iy);
          const imaginaryPositionZ = positionsJ[2].plus(iz);

          const rij = [
            positionsI[0].minus(imaginaryPositionX),
            positionsI[1].minus(imaginaryPositionY),
            positionsI[2].minus(imaginaryPositionZ),
          ];

          const bij = rij[0].times(uij[0]).plus(rij[1].times(uij[1])).plus(rij[2].times(uij[2]));
          const cij = rij[0].times(rij[0]).plus(rij[1].times(rij[1])).plus(rij[2].times(rij[2])).minus(sigmaSquared);

          if (bij.lt(0)) {
            const discriminant = bij.times(bij).minus(uij2.times(cij));
            if (discriminant.gt(0)) {
              const time = bij.neg().minus(discriminant.sqrt()).div(uij2);
              if (time.lt(this.collisionTime[i][j])) {
                this.collisionTime[i][j] = time;
              }
            }
          }
        }
      }
    }
  }

  writeResults() {
    const fs = require('fs');
    const XLSX = require('xlsx');
    const wb = XLSX.utils.book_new();
    const ws = XLSX.utils.json_to_sheet(this.results);
    XLSX.utils.book_append_sheet(wb, ws, 'Table Data');
    fs.writeFileSync(`./output/results_${this.nSpheres.toNumber()}_${this.rVolume.toFixed(2)}_${this.nCollisions}.xlsx`, XLSX.write(wb, {type: 'buffer'}));
    console.log('Table data written to file');
  }

}
