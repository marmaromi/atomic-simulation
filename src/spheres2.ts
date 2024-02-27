import {spheresResultsPropertiesModel} from './models';
import {Calc} from './Calc';

export class Spheres2 {
  nSpheres: number;
  rVolume: number;
  nCollisions: number;
  sigma: number;
  positions: number[][] = [];
  velocities: number[][] = [];
  collisionTime: number[][] = [];
  tCol: number;
  iCol: number;
  jCol: number;
  deltaVel: number[];
  momentum: number[];
  kineticEnergy: number;
  vSample: number;
  orderP: number;
  dv: number;
  pv0: number;
  totalDV = 0;
  totalTCOL = 0;
  results: spheresResultsPropertiesModel[] = [];

// Constants
  boxLength: number = 1;


  constructor(nSpheres: number, rVolume: number, nCollisions: number) {
    this.setInput(nSpheres, rVolume, nCollisions);
    this.validateInput();
    this.computeDiameter();
    this.assignPositions();
    this.assignVelocitiesZeroMom();
    this.writeInitial();
    this.updateCollisionsTable('init');

    //start collision loop
    for (let i = 0; i < this.nCollisions; i++) {
      this.retrieveCollisionsInfo();
      this.advanceSimulation();
      this.computeProperties();
      this.writeProperties(i + 1);
      this.updateCollisionsTable('update');

      if (i > 5000) { // skip first n collisions
        this.totalDV += this.dv;
        this.totalTCOL += this.tCol;
      }
    }
    this.pv0 = this.rVolume * (this.totalDV * (this.sigma / (this.totalTCOL * this.nSpheres * 3)) + 1);
    this.writeResults();
  }

  setInput (nSpheres: number, rVolume: number, nCollisions: number) {
    this.nSpheres = nSpheres;
    this.rVolume = rVolume;
    this.nCollisions = nCollisions;
  }

  validateInput() {
    const nInt = Math.pow(this.nSpheres / 4, 1 / 3);
    const nIntRounded = Math.round(nInt * 1e12) / 1e12;
    if (!Number.isInteger(nIntRounded)) {
      throw Error(`Wrong nSpheres!, nSpheres must be a perfect cube of 4. nSpheres = ${this.nSpheres}`);
    }
    if (this.rVolume < 1) {
      throw Error(`Wrong rVolume!, rVolume must be >= 1. rVolume = ${this.rVolume}`);
    }
  }

  computeDiameter() {
    this.sigma = Math.pow(Math.sqrt(2) / (this.nSpheres * this.rVolume), 1/3) * this.boxLength;
  }

  assignPositions() {
    const nInt = Math.pow(this.nSpheres / 4, 1 / 3);
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
    this.positions = positions;
  }

  assignVelocitiesZeroMom() {
    const speed = Math.sqrt(3);
    const twoPi = 2 * Math.PI;
    const velocities = [];
    let vSumX = 0;
    let vSumY = 0;
    let vSumZ = 0;
    for (let i = 0; i < this.nSpheres; i++) {
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
    const vAvgX = vSumX / this.nSpheres;
    const vAvgY = vSumY / this.nSpheres;
    const vAvgZ = vSumZ / this.nSpheres;

    for (let i = 0; i < this.nSpheres; i++) {
      velocities[i][0] -= vAvgX;
      velocities[i][1] -= vAvgY;
      velocities[i][2] -= vAvgZ;
    }
    this.velocities = velocities;
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
    let tCol = Infinity;
    let iCol = -1;
    let jCol = -1;

    for (let i = 0; i < n - 1; i++) {
      for (let j = i + 1; j < n; j++) {
        if (this.collisionTime[i][j] < tCol) {
          tCol = this.collisionTime[i][j];
          iCol = i;
          jCol = j;
        }
      }
    }

    this.tCol = tCol;
    this.iCol = iCol;
    this.jCol = jCol;
  }

  advanceSimulation() {
    for (let i = 0; i < this.nSpheres; i++) {
      //update positions
      this.positions[i] = [
        this.positions[i][0] + this.velocities[i][0] * this.tCol,
        this.positions[i][1] + this.velocities[i][1] * this.tCol,
        this.positions[i][2] + this.velocities[i][2] * this.tCol,
      ];

      // If the new position is outside the simulation wall
      this.positions[i][0] = this.positions[i][0] % 1;
      this.positions[i][1] = this.positions[i][1] % 1;
      this.positions[i][2] = this.positions[i][2] % 1;
    }

    //update velocities of next collision
    const uij = [
      this.velocities[this.iCol][0] - this.velocities[this.jCol][0],
      this.velocities[this.iCol][1] - this.velocities[this.jCol][1],
      this.velocities[this.iCol][2] - this.velocities[this.jCol][2],
    ];
    let minDistance = Infinity;
    let bij;
    let rij: number[] = [];

    for (let ix = -1; ix <= 1; ix++) {
      for (let iy = -1; iy <= 1; iy++) {
        for (let iz = -1; iz <= 1; iz++) {
          const translate = [ix, iy, iz];
          const imaginaryPosition = [
            this.positions[this.jCol][0] + translate[0],
            this.positions[this.jCol][1] + translate[1],
            this.positions[this.jCol][2] + translate[2],
          ];
          const rijImaginary = [
            this.positions[this.iCol][0] - imaginaryPosition[0],
            this.positions[this.iCol][1] - imaginaryPosition[1],
            this.positions[this.iCol][2] - imaginaryPosition[2],
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
    const sigmaSquared = Math.pow(this.sigma, 2);
    this.deltaVel = [
      -(bij / sigmaSquared) * rij[0],
      -(bij / sigmaSquared) * rij[1],
      -(bij / sigmaSquared) * rij[2]
    ];
    this.velocities[this.iCol] = [
      this.velocities[this.iCol][0] + this.deltaVel[0],
      this.velocities[this.iCol][1] + this.deltaVel[1],
      this.velocities[this.iCol][2] + this.deltaVel[2],
    ];
    this.velocities[this.jCol] = [
      this.velocities[this.jCol][0] - this.deltaVel[0],
      this.velocities[this.jCol][1] - this.deltaVel[1],
      this.velocities[this.jCol][2] - this.deltaVel[2],
    ];
  }

  computeProperties() {
    this.momentum = [0, 0, 0];
    this.kineticEnergy = 0;
    this.vSample = 0;
    this.orderP = 0;

    const a = this.boxLength / Math.pow(this.nSpheres / 4, 1 / 3);

    for (let i = 0; i < this.nSpheres; i++) {
      this.momentum = [
        this.momentum[0] + this.velocities[i][0],
        this.momentum[1] + this.velocities[i][1],
        this.momentum[2] + this.velocities[i][2],
      ];
      const v = Math.sqrt(Calc.dotProduct2(this.velocities[i], this.velocities[i]));
      this.kineticEnergy += 0.5 * v ** 2;

      if (1 <= v && v < 2) {
        this.vSample += 1;
      }

      this.orderP += Math.cos(4 * Math.PI * this.positions[i][0] / a) +
        Math.cos(4 * Math.PI * this.positions[i][1] / a) +
        Math.cos(4 * Math.PI * this.positions[i][2] / a);
    }

    this.momentum = [this.momentum[0] / this.nSpheres, this.momentum[1] / this.nSpheres, this.momentum[2] / this.nSpheres];
    this.kineticEnergy /= this.nSpheres;
    this.orderP /= 3 * this.nSpheres;
    this.dv = Math.sqrt(Calc.dotProduct2(this.deltaVel, this.deltaVel));
  }

  writeProperties(nthCollision: number) {
    this.results.push({
      n: nthCollision,
      momX: this.momentum[0].toFixed(3),
      momY: this.momentum[1].toFixed(3),
      momZ: this.momentum[2].toFixed(3),
      kinE: this.kineticEnergy.toFixed(3),
      tCol: this.tCol,
      vSample: this.vSample,
      orderP: this.orderP.toFixed(3),
      dv: this.dv,
    });
  }

  updateCollisionsTable(role: string) {
    if (role === 'init') {
      this.collisionTime = new Array(this.nSpheres);
      for (let index = 0; index < this.nSpheres; index++) {
        this.collisionTime[index] = new Array(this.nSpheres).fill(Infinity);
      }
    }

    for (let i = 1; i <= this.nSpheres - 1; i++) {
      for (let j = i + 1; j <= this.nSpheres; j++) {
        const k = i - 1;
        const l = j - 1;
        if (role === 'init') {
          this.loopCollisionsTable(k, l);
        } else if (role === 'update') {
          this.collisionTime[k][l] -= this.tCol;
          if (k === this.iCol || k === this.jCol || l === this.iCol || l === this.jCol) {
            this.collisionTime[k][l] = Infinity;
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
    const sigmaSquared = this.sigma * this.sigma;
    const collisionTimeIJ = this.collisionTime[i][j];

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
    fs.writeFileSync(`./output/results_${this.nSpheres}_${this.rVolume.toFixed(2)}_${this.nCollisions}.xlsx`, XLSX.write(wb, {type: 'buffer'}));
    console.log('Table data written to file');
  }

}

