import {Calc} from './Calc';
import {spheresResultsPropertiesModel} from './models';

export class SpheresCopilot {

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
  totalDV: number = 0;
  totalTCOL: number = 0;
  results: spheresResultsPropertiesModel[] = [];

  boxLength = 1;
  PI = Math.acos(-1);

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

  setInput(nSpheres: number, rVolume: number, nCollisions: number) {
    this.nSpheres = nSpheres;
    this.rVolume = rVolume;
    this.nCollisions = nCollisions;
  }

  validateInput() {
    const nInt: number = Math.pow(this.nSpheres / 4, 1 / 3);
    const nIntRounded = Math.round(nInt * 1e12) / 1e12;
    console.log(nIntRounded);
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

    for (let i = 0; i < nInt; i++) {
      for (let j = 0; j < nInt; j++) {
        for (let k = 0; k < nInt; k++) {
          const translation = [a * i, a * j, a * k];
          this.positions.push(translation);

          const trans1 = [a / 2, a / 2, 0];
          this.positions.push(translation.map((loc, index) => loc + trans1[index]));

          const trans2 = [a / 2, 0, a / 2];
          this.positions.push(translation.map((loc, index) => loc + trans2[index]));

          const trans3 = [0, a / 2, a / 2];
          this.positions.push(translation.map((loc, index) => loc + trans3[index]));
        }
      }
    }
  }

  assignVelocitiesZeroMom() {
    const speed = Math.sqrt(3);
    for (let i = 1; i <= this.nSpheres; i++) {
      const u = Math.random();
      const v = Math.random();
      const theta = this.PI * 2 * u;
      const phi = Math.acos(2 * v - 1);
      this.velocities.push([
        speed * Math.sin(phi) * Math.cos(theta),
        speed * Math.sin(phi) * Math.sin(theta),
        speed * Math.cos(phi),
      ]);
    }
    const vSum = this.velocities.reduce((a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
    const vAvg = vSum.map(i => i / this.nSpheres);
    this.velocities = this.velocities.map(vel => [vel[0] - vAvg[0], vel[1] - vAvg[1], vel[2] - vAvg[2]]);
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
    this.tCol = Infinity;

    for (let i = 1; i <= n - 1; i++) {
      for (let j = i + 1; j <= n; j++) {
        if (this.collisionTime[i - 1][j - 1] < this.tCol) {
          this.tCol = this.collisionTime[i - 1][j - 1];
          this.iCol = i - 1;
          this.jCol = j - 1;
        }
      }
    }
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
    let uij = [
      this.velocities[this.iCol][0] - this.velocities[this.jCol][0],
      this.velocities[this.iCol][1] - this.velocities[this.jCol][1],
      this.velocities[this.iCol][2] - this.velocities[this.jCol][2],
    ];
    let minDistance = Infinity;
    let bij: number;
    let rij: number[];

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
    this.deltaVel = [-(bij / this.sigma ** 2) * rij[0], -(bij / this.sigma ** 2) * rij[1], -(bij / this.sigma ** 2) * rij[2]];
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

      if (v > 1 && v < 2) {
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
      for (let index = 0; index < this.collisionTime.length; index++) {
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
          this.collisionTime[k][l] = this.collisionTime[k][l] - this.tCol;
          if (k === this.iCol || k === this.jCol || l === this.iCol || l === this.jCol) {
            this.collisionTime[k][l] = Infinity;
            this.loopCollisionsTable(k, l);
          }
        }
      }
    }
  }

  loopCollisionsTable(i: number, j: number) {
    const uij = [
      this.velocities[i][0] - this.velocities[j][0],
      this.velocities[i][1] - this.velocities[j][1],
      this.velocities[i][2] - this.velocities[j][2],
    ];
    for (let ix = -1; ix <= 1; ix++) {
      for (let iy = -1; iy <= 1; iy++) {
        for (let iz = -1; iz <= 1; iz++) {
          const translate = [ix, iy, iz];
          const imaginaryPosition = [
            this.positions[j][0] + translate[0],
            this.positions[j][1] + translate[1],
            this.positions[j][2] + translate[2],
          ];
          const rij = [
            this.positions[i][0] - imaginaryPosition[0],
            this.positions[i][1] - imaginaryPosition[1],
            this.positions[i][2] - imaginaryPosition[2],
          ];
          const bij = Calc.dotProduct2(rij, uij);
          const cij = Calc.dotProduct2(rij, rij) - this.sigma ** 2;

          if (bij < 0) {
            const uij2 = Calc.dotProduct2(uij, uij);
            const discriminant = bij * bij - uij2 * cij;
            if (discriminant > 0) {
              const time = -bij - Math.sqrt(discriminant) / uij2;
              if (time < this.collisionTime[i][j]) {
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
