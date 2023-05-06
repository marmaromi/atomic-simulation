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

  boxLength = 1;
  PI = Decimal.acos(-1);

  constructor(nSpheres: number, rVolume: number, nCollisions: number) {
    this.setInput(nSpheres, rVolume, nCollisions);
    this.validateInput();
    this.computeDiameter();
    this.assignPositions();
    this.assignVelocitiesZeroMom();
    this.writeInitial();
    this.initializeCollisionsTable();
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
    const nInt: Decimal = this.nSpheres.dividedBy(4).pow(1 / 3);
    if (!Number.isInteger(nInt.toNumber())) {
      throw Error('Wrong nSpheres!');
    }
    if (this.rVolume < 1) {
      throw Error('Wrong rVolume!');
    }
    // console.groupEnd();
  }

  computeDiameter() {
    console.group('computeDiameter');
    this.sigma = Decimal.sqrt(2).dividedBy(this.nSpheres.mul(this.rVolume)).pow(1 / 3).mul(this.boxLength);
    // console.log('Decimal sigma', this.sigma);
    // console.log('math sigma',Math.pow(Math.sqrt(2) / (this.nSpheres.toNumber() * this.rVolume), 1 / 3) * this.boxLength);
    console.groupEnd();
  }

  assignPositions() {
    // console.group('assignPositions');

    // const nInt = Math.pow(this.nSpheres / 4, 1 / 3);
    // const a = 1 / nInt;
    const nInt = this.nSpheres.dividedBy(4).pow(1 / 3);
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
    // let vSum = []
    for (let i = 1; i <= this.nSpheres.toNumber(); i++) {
      // const u = Math.random();
      // const v = Math.random();
      // const theta = 2 * Math.PI * u;
      // const phi = Math.acos(2 * v - 1);
      // this.velocities.push([
      //   speed * Math.sin(phi) * Math.cos(theta),
      //   speed * Math.sin(phi) * Math.sin(theta),
      //   speed * Math.cos(phi),
      // ]);
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
    // const vSum = this.velocities.reduce((a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
    const vSum = this.velocities.reduce((a, b) => [a[0].add(b[0]), a[1].add(b[1]), a[2].add(b[2])]);
    // console.log('vSum: ', vSum);

    // const vAvg = vSum.map(i => i / this.nSpheres);
    const vAvg = vSum.map(i => i.dividedBy(this.nSpheres));
    // console.log('vAvg', vAvg);

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
    console.log('Initial positions of the spheres: ', this.positions);
    console.log('Initial velocities of the spheres: ', this.velocities);
    console.groupEnd();
  }

  initializeCollisionsTable() {
    console.group('initializeCollisionsTable: ');
    this.collisionTime = new Array(this.nSpheres.toNumber());
    for (let index = 0; index < this.collisionTime.length; index++) {
      this.collisionTime[index] = new Array(this.nSpheres.toNumber()).fill(new Decimal(Infinity));
    }
    for (let i = 0; i < this.nSpheres.toNumber() - 1; i++) {
      for (let j = i + 1; j < this.nSpheres.toNumber(); j++) {
        // console.log(i,j);
        const uij = [this.velocities[i][0].minus(this.velocities[j][0]), this.velocities[i][1].minus(this.velocities[j][1]), this.velocities[i][2].minus(this.velocities[j][2])];
        // console.log('uij',uij);
        for (let ix = -1; ix <= 1; ix++) {
          for (let iy = -1; iy <= 1; iy++) {
            for (let iz = -1; iz <= 1; iz++) {
              const translate = [ix, iy, iz];
              const imaginaryPosition = [this.positions[j][0].add(translate[0]), this.positions[j][1].add(translate[1]), this.positions[j][2].add(translate[2])];
              // console.log('translate',translate);
              // console.log('imaginaryPosition',imaginaryPosition);
              const rij = [this.positions[i][0].minus(imaginaryPosition[0]), this.positions[i][1].minus(imaginaryPosition[1]), this.positions[i][2].minus(imaginaryPosition[2])];
              // console.log('rij',rij);
              const bij = Calc.dotProduct(rij, uij);
              // console.log(rij,uij,bij);
              // const cij = Calc.dotProduct(rij, rij) - Math.pow(this.sigma, 2);
              const cij = Calc.dotProduct(rij, rij).minus(this.sigma.pow(2));
              if (bij < new Decimal(0)) {
                const uij2 = Calc.dotProduct(uij, uij);
                const disc = bij.mul(bij).minus(uij2.mul(cij));
                if (disc > new Decimal(0)) {
                  // console.log('uij: ',uij);
                  // console.log('uij','bij');
                  // console.log(uij,bij);
                  // console.log('bij', bij);
                  // console.log('disc', disc);
                  // console.log('uij2: ',uij2);
                  // const time = (-bij - Math.sqrt(disc)) / uij2;
                  const time = new Decimal(-bij.minus(Decimal.sqrt(disc).dividedBy(uij2)));
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
        // console.group('time for collision');
        // console.log('[i,j]:               ',[i,j]);
        // console.log('collisionTime[i][j]: ', this.collisionTime[i][j]);
        console.groupEnd()
      }
    }
    console.log(this.collisionTime);
    console.groupEnd();
  }


}