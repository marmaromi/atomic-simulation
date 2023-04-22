export class Spheres {

  nSpheres: number;
  rVolume: number;
  nCollisions: number;
  sigma: number;
  positions: number[][] = [];
  velocities: number[][] = [];
  collisionTime: number[][] = [];

  boxLength = 1;


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
    this.nSpheres = nSpheres;
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
    const nInt = Math.pow(this.nSpheres / 4, 1 / 3);
    // console.log('nInt: ', nInt);
    if (!Number.isInteger(nInt)) {
      throw Error('Wrong nSpheres!');
    }
    if (this.rVolume < 1) {
      throw Error('Wrong rVolume!');
    }
    // console.groupEnd();
  }

  computeDiameter() {
    // console.group('computeDiameter');
    this.sigma = Math.pow(Math.sqrt(2) / (this.nSpheres / this.rVolume), 1 / 3) * this.boxLength;
    // console.log('sigma', this.sigma);
    // console.groupEnd();
  }

  assignPositions() {
    // console.group('assignPositions');

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
    // console.log(this.positions);
    // console.groupEnd();
  }

  assignVelocitiesZeroMom() {
    // console.group('assignVelocitiesZeroMom');
    const speed = Math.sqrt(3);
    // let vSum = []
    for (let i = 1; i <= this.nSpheres; i++) {
      const u = Math.random();
      const v = Math.random();
      const theta = 2 * Math.PI * u;
      const phi = Math.acos(2 * v - 1);
      this.velocities.push([
        speed * Math.sin(phi) * Math.cos(theta),
        speed * Math.sin(phi) * Math.sin(theta),
        speed * Math.cos(phi)]);
    }
    // console.log('velocities: ', this.velocities);
    const vSum = this.velocities.reduce((a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
    // console.log('vSum: ', vSum);

    const vAvg = vSum.map(i => i / this.nSpheres);
    // console.log('vAvg', vAvg);

    this.velocities = this.velocities.map(arr => [arr[0] - vAvg[0], arr[1] - vAvg[1], arr[2] - vAvg[2]]);
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
    this.collisionTime = new Array(this.nSpheres);
    for (let index = 0; index < this.collisionTime.length; index++) {
      this.collisionTime[index] = new Array(this.nSpheres - 1).fill(Number.MAX_VALUE);
    }
    for (let i = 0; i < this.nSpheres; i++) {
      for (let j = i + 1; j < this.nSpheres - 1; j++) {
        const uij = [this.velocities[i][0] - this.velocities[j][0], this.velocities[i][1] - this.velocities[j][1], this.velocities[i][2] - this.velocities[j][2]];
        // console.log(uij);
        for (let ix = -1; ix <= 1; ix++) {
          for (let iy = -1; iy <= 1; iy++) {
            for (let iz = -1; iz <= 1; iz++) {
              const translate = [ix, iy, iz];
              const imaginaryPosition = [this.positions[j][0] - translate[0], this.positions[j][1] - translate[1], this.positions[j][2] - translate[2]];
              const rij = [this.positions[i][0] - imaginaryPosition[0], this.positions[i][1] - imaginaryPosition[1], this.positions[i][2] - imaginaryPosition[2]];
              const bij = this.dotProduct(rij, uij);
              const cij = this.dotProduct(rij, rij) - Math.pow(this.sigma, 2);
              if (bij < 0) {
                const uij2 = this.dotProduct(uij, uij);
                const disc = bij * bij - uij2 * cij;
                if (disc > 0) {
                  const time = (-bij - Math.sqrt(disc)) / uij2;
                  if (time < this.collisionTime[i][j]) {
                    this.collisionTime[i][j] = time;
                  }
                }
              }
            }
          }
        }
      }
    }
    console.log(this.collisionTime);
    console.groupEnd();
  }


  dotProduct(a: number[], b: number[]): number {
    return a.map((x, i) => a[i] * b[i]).reduce((m, n) => m + n);
  }

}