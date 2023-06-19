import {Decimal} from 'decimal.js';

export class Calc {
  static dotProduct(a: Decimal[], b: Decimal[]): Decimal {
    return a.map((x, i) => a[i].mul(b[i])).reduce((m, n) => m.add(n));
  }
}