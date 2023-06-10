import {Calc} from './Calc';
import {Decimal} from 'decimal.js';

export class Spheres {

  nSpheres: Decimal;
  rVolume: number;
  nCollisions: number;
  sigma: Decimal;
  positions: Decimal[][] = [];
  velocities: Decimal[][] = [];
  collisionTime: Decimal[][] = [];
  tCol: Decimal;
  iCol: number;
  jCol: number;
  deltaVel: Decimal[];

  momentum: Decimal[];
  kineticEnergy: Decimal;
  vSample: Decimal;
  orderP: Decimal;
  dv: Decimal;

  results: any = [];


  boxLength = new Decimal(1);
  PI = Decimal.acos(-1);

  constructor(nSpheres: number, rVolume: number, nCollisions: number) {
    this.setInput(nSpheres, rVolume, nCollisions);
    this.validateInput();
    this.computeDiameter();
    this.assignPositions();
    this.assignVelocitiesZeroMom();
    this.writeInitial();
    this.initializeCollisionsTable();

    //start collision loop
    for (let i = 0; i < this.nCollisions; i++) {
      this.retrieveCollisionsInfo();
      this.advanceSimulation();
      this.computeProperties();
      this.writeProperties(i + 1);
      this.updateCollisionsTable();
    }
    this.writeResults();
  }


  setInput(nSpheres: number, rVolume: number, nCollisions: number) {
    // console.group('setInput');
    this.nSpheres = new Decimal(nSpheres);
    this.rVolume = rVolume;
    this.nCollisions = nCollisions;
    // console.log('nSpheres: ', this.nSpheres);
    // console.log('rVolume: ', this.rVolume);
    // console.log('nCollisions: ', this.nCollisions);
    //
    // console.groupEnd();
  }

  validateInput() {
    // console.group('validateInput');
    const nInt: Decimal = this.nSpheres.dividedBy(4).pow(Decimal.div(1, 3));
    if (!Number.isInteger(nInt.toNumber())) {
      throw Error('Wrong nSpheres!');
    }
    if (this.rVolume < 1) {
      throw Error('Wrong rVolume!');
    }
    // console.groupEnd();
  }

  computeDiameter() {
    // console.group('computeDiameter');
    this.sigma = Decimal.sqrt(2).dividedBy(this.nSpheres.mul(this.rVolume)).pow(Decimal.div(1, 3)).mul(this.boxLength);
    // console.log('sigma', this.sigma);
    // console.groupEnd();
  }

  assignPositions() {
    // console.group('assignPositions');
    const nInt = this.nSpheres.dividedBy(4).pow(Decimal.div(1, 3));
    const a = new Decimal(1).dividedBy(nInt);

    for (let i = 0; i < nInt.toNumber(); i++) {
      for (let j = 0; j < nInt.toNumber(); j++) {
        for (let k = 0; k < nInt.toNumber(); k++) {
          const translation = [a.mul(i), a.mul(j), a.mul(k)];
          this.positions.push(translation);

          const trans1 = [a.dividedBy(2), a.dividedBy(2), 0];
          this.positions.push(translation.map((loc, index) => loc.add(trans1[index])));

          const trans2 = [a.dividedBy(2), 0, a.dividedBy(2)];
          this.positions.push(translation.map((loc, index) => loc.add(trans2[index])));

          const trans3 = [0, a.dividedBy(2), a.dividedBy(2)];
          this.positions.push(translation.map((loc, index) => loc.add(trans3[index])));
        }
      }
    }
    // console.log(this.positions);
    // console.groupEnd();
  }

  assignVelocitiesZeroMom() {
    // this.velocities = [
    //   [new Decimal(-0.64560687459701493), new Decimal(-0.88029736677446080), new Decimal(-0.54225736800773905)],
    //   [new Decimal(-0.59679349070925591), new Decimal(-0.093737658514908229), new Decimal(-1.1491893834322298)],
    //   [new Decimal(1.4577059995277539), new Decimal(-0.32175617027540171), new Decimal(-0.18060048926727768)],
    //   [new Decimal(-0.21530563422148283), new Decimal(1.2957911955647707), new Decimal(1.8720472407072464)],
    // ];
    // return;
    // console.group('assignVelocitiesZeroMom');
    // const speed = Math.sqrt(3);
    const speed = Decimal.sqrt(3);
    for (let i = 1; i <= this.nSpheres.toNumber(); i++) {
      const u = Decimal.random();
      const v = Decimal.random();
      const theta = this.PI.mul(2).mul(u);
      const phi = Decimal.acos(v.mul(2).minus(1));
      this.velocities.push([
        speed.mul(Decimal.sin(phi)).mul(Decimal.cos(theta)),
        speed.mul(Decimal.sin(phi)).mul(Decimal.sin(theta)),
        speed.mul(Decimal.cos(phi)),
      ]);
    }
    // console.log('velocities: ', this.velocities);
    const vSum = this.velocities.reduce((a, b) => [a[0].add(b[0]), a[1].add(b[1]), a[2].add(b[2])]);
    // console.log('vSum: ', vSum);

    // const vAvg = vSum.map(i => i / this.nSpheres);
    const vAvg = vSum.map(i => i.div(this.nSpheres));

    // this.velocities = this.velocities.map(arr => [arr[0] - vAvg[0], arr[1] - vAvg[1], arr[2] - vAvg[2]]);
    this.velocities = this.velocities.map(arr => [arr[0].minus(vAvg[0]), arr[1].minus(vAvg[1]), arr[2].minus(vAvg[2])]);
    // console.log('velocities zero momentum: ', this.velocities);

    // console.groupEnd();

  }

  writeInitial() {
    console.group('Input received:');
    console.log('Number of spheres: ', this.nSpheres);
    console.log('Reduced volume: ', this.rVolume);
    console.log('Number of collisions: ', this.nCollisions);
    // console.log('Initial positions of the spheres: ', this.positions);
    // console.log('Initial velocities of the spheres: ', this.velocities);
    console.groupEnd();
    console.log('----------------------------------------');
  }

  initializeCollisionsTable() {
    // console.group('initializeCollisionsTable: ');
    this.collisionTime = new Array(this.nSpheres.toNumber());
    for (let index = 0; index < this.collisionTime.length; index++) {
      this.collisionTime[index] = new Array(this.nSpheres.toNumber()).fill(new Decimal(Infinity));
    }
    for (let i = 0; i < this.nSpheres.toNumber() - 1; i++) {
      for (let j = i + 1; j < this.nSpheres.toNumber(); j++) {
        // console.log(i,j);
        this.loopCollisionsTable(i, j);
        // console.group('time for collision');
        // console.log('[i,j]:               ',[i,j]);
        // console.log('collisionTime[i][j]: ', this.collisionTime[i][j]);
      }
    }
    // console.log(this.collisionTime);
    // console.groupEnd()
  }

  retrieveCollisionsInfo() {
    // console.group('retrieveCollisionsInfo: ')
    const n = this.collisionTime.length;
    this.tCol = new Decimal(Infinity);

    for (let i = 1; i <= n - 1; i++) {
      for (let j = i + 1; j <= n; j++) {
        if (this.collisionTime[i-1][j-1] < this.tCol) {
          this.tCol = this.collisionTime[i-1][j-1];
          this.iCol = i-1;
          this.jCol = j-1;
        }
      }
    }
    // console.log('tCol: ', this.tCol);
    // console.log('iCol: ', this.iCol);
    // console.log('jCol: ', this.jCol);
    // console.groupEnd();
  }

  advanceSimulation() {
    for (let i = 0; i < this.nSpheres.toNumber(); i++) {
      //update positions
      this.positions[i] = [
        this.positions[i][0].add(this.velocities[i][0].times(this.tCol)),
        this.positions[i][1].add(this.velocities[i][1].times(this.tCol)),
        this.positions[i][2].add(this.velocities[i][2].times(this.tCol)),
      ];

      // If the new position is outside the simulation wall
      this.positions[i][0] = this.positions[i][0].mod(1);
      this.positions[i][1] = this.positions[i][1].mod(1);
      this.positions[i][2] = this.positions[i][2].mod(1);
    }

    //update velocities of next collision
    let uij = [
      this.velocities[this.iCol][0].minus(this.velocities[this.jCol][0]),
      this.velocities[this.iCol][1].minus(this.velocities[this.jCol][1]),
      this.velocities[this.iCol][2].minus(this.velocities[this.jCol][2]),
    ];
    let minDistance = new Decimal(Infinity);
    let bij: Decimal;
    let rij: Decimal[];

    for (let ix = -1; ix <= 1; ix++) {
      for (let iy = -1; iy <= 1; iy++) {
        for (let iz = -1; iz <= 1; iz++) {
          const translate = [ix, iy, iz];
          const imaginaryPosition = [
            this.positions[this.jCol][0].add(translate[0]),
            this.positions[this.jCol][1].add(translate[1]),
            this.positions[this.jCol][2].add(translate[2]),
          ];
          const rijImaginary = [
            this.positions[this.iCol][0].minus(imaginaryPosition[0]),
            this.positions[this.iCol][1].minus(imaginaryPosition[1]),
            this.positions[this.iCol][2].minus(imaginaryPosition[2]),
          ];
          const distance = Calc.dotProduct(rijImaginary, rijImaginary);
          if (distance < minDistance) {
            rij = rijImaginary;
            bij = Calc.dotProduct(rijImaginary, uij);
            minDistance = distance;
          }
        }
      }
    }
    this.deltaVel = [(bij.div(this.sigma.pow(2))).mul(rij[0]).neg(), (bij.div(this.sigma.pow(2))).mul(rij[1]).neg(), (bij.div(this.sigma.pow(2))).mul(rij[2]).neg()];
    this.velocities[this.iCol] = [
      this.velocities[this.iCol][0].add(this.deltaVel[0]),
      this.velocities[this.iCol][1].add(this.deltaVel[1]),
      this.velocities[this.iCol][2].add(this.deltaVel[2]),
    ];
    this.velocities[this.jCol] = [
      this.velocities[this.jCol][0].minus(this.deltaVel[0]),
      this.velocities[this.jCol][1].minus(this.deltaVel[1]),
      this.velocities[this.jCol][2].minus(this.deltaVel[2]),
    ];

    // console.log('deltaVel', this.deltaVel);
    // console.log('velocities[iCol]', this.velocities[this.iCol]);
    // console.log('velocities[jCol]', this.velocities[this.jCol]);

  }

  computeProperties() {
    this.momentum = [new Decimal(0), new Decimal(0), new Decimal(0)];
    this.kineticEnergy = new Decimal(0);
    this.vSample = new Decimal(0);
    this.orderP = new Decimal(0);

    const a = this.boxLength.div((this.nSpheres.div(4)).pow(Decimal.div(1, 3)));

    for (let i = 0; i < this.nSpheres.toNumber(); i++) {
      this.momentum = [
        this.momentum[0].add(this.velocities[i][0]),
        this.momentum[1].add(this.velocities[i][1]),
        this.momentum[2].add(this.velocities[i][2]),
      ];
      const v = Decimal.sqrt(Calc.dotProduct(this.velocities[i], this.velocities[i]));
      this.kineticEnergy = this.kineticEnergy.add(Decimal.mul(0.5, v.pow(2)));

      if (v.greaterThan(1) && v.lessThan(2)) {
        this.vSample = this.vSample.add(1);
      }

      this.orderP = this.orderP
        .add(Decimal.cos(this.PI.mul(4).mul(this.positions[i][0]).div(a)))
        .add(Decimal.cos(this.PI.mul(4).mul(this.positions[i][1]).div(a)))
        .add(Decimal.cos(this.PI.mul(4).mul(this.positions[i][2]).div(a)));
    }

    this.momentum = [this.momentum[0].div(this.nSpheres), this.momentum[1].div(this.nSpheres), this.momentum[2].div(this.nSpheres)];
    this.kineticEnergy = this.kineticEnergy.div(this.nSpheres);
    this.orderP = this.orderP.div(this.nSpheres.mul(3));
    this.dv = Decimal.sqrt(Calc.dotProduct(this.deltaVel, this.deltaVel));

    // console.group('computeProperties');
    // console.log('momentum', this.momentum);
    // console.log('kineticEnergy', this.kineticEnergy);
    // console.log('vSample', this.vSample);
    // console.log('orderP', this.orderP);
    // console.log('dv', this.dv);
    // console.groupEnd();
  }

  writeProperties(nthCollision: number) {
    this.results.push({
      n: nthCollision,
      mom0: this.momentum[0].toFixed(3),
      mom1: this.momentum[1].toFixed(3),
      mom2: this.momentum[2].toFixed(3),
      kin: this.kineticEnergy.toFixed(3),
      tCol: this.tCol.toNumber(),
      vSample: this.vSample.toNumber(),
      orderP: this.orderP.toFixed(3),
      dv: this.dv.toNumber()});
  }

  updateCollisionsTable() {
    for (let i = 0; i < this.nSpheres.toNumber() - 1; i++) {
      for (let j = i + 1; j < this.nSpheres.toNumber(); j++) {
        this.collisionTime[i][j] = this.collisionTime[i][j].minus(this.tCol);
        if (i === this.iCol || i === this.jCol || j === this.iCol || j === this.jCol) {
          this.collisionTime[i][j] = new Decimal(Infinity);
          this.loopCollisionsTable(i, j);
        }
      }
    }
  }

  writeResults() {
    console.table(this.results);
  }

  loopCollisionsTable(i: number, j: number) {
    const uij = [
      this.velocities[i][0].minus(this.velocities[j][0]),
      this.velocities[i][1].minus(this.velocities[j][1]),
      this.velocities[i][2].minus(this.velocities[j][2]),
    ];
    for (let ix = -1; ix <= 1; ix++) {
      for (let iy = -1; iy <= 1; iy++) {
        for (let iz = -1; iz <= 1; iz++) {
          const translate = [ix, iy, iz];
          const imaginaryPosition = [
            this.positions[j][0].add(translate[0]),
            this.positions[j][1].add(translate[1]),
            this.positions[j][2].add(translate[2]),
          ];
          // console.log('translate',translate);
          // console.log('imaginaryPosition',imaginaryPosition);
          const rij = [
            this.positions[i][0].minus(imaginaryPosition[0]),
            this.positions[i][1].minus(imaginaryPosition[1]),
            this.positions[i][2].minus(imaginaryPosition[2]),
          ];
          // console.log('rij',rij);
          const bij = Calc.dotProduct(rij, uij);
          // console.log(rij,uij,bij);
          // const cij = Calc.dotProduct(rij, rij) - Math.pow(this.sigma, 2);
          const cij = Calc.dotProduct(rij, rij).minus(this.sigma.pow(2));
          if (bij.lessThan(new Decimal(0))) {
            const uij2 = Calc.dotProduct(uij, uij);
            const discriminant = bij.mul(bij).minus(uij2.mul(cij));
            if (discriminant.greaterThan(new Decimal(0))) {
              // console.log('uij: ',uij);
              // console.log('uij','bij');
              // console.log(uij,bij);
              // console.log('bij', bij);
              // console.log('discriminant', discriminant);
              // console.log('uij2: ',uij2);
              // const time = (-bij - Math.sqrt(discriminant)) / uij2;
              const time = new Decimal(-bij.minus(Decimal.sqrt(discriminant).dividedBy(uij2)));
              // console.group('time for collision');
              // console.log('time:                ', time);
              // console.log('collisionTime[i][j]: ', this.collisionTime[i][j]);
              // console.groupEnd()
              if (time < this.collisionTime[i][j]) {
                // console.group('time for collision');
                // console.log('time:                ', time);
                // console.log('collisionTime[i][j]: ', this.collisionTime[i][j]);
                // console.groupEnd()
                this.collisionTime[i][j] = time;
              }
              // console.log('this.collisionTime[i][j]',this.collisionTime[i][j]);
            }
          }
        }
      }
    }
  }
}
