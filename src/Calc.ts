import {Decimal} from 'decimal.js';
import {Big} from 'big.js';


export class Calc {
  static dotProduct(a: Decimal[], b: Decimal[]): Decimal {
    return a.map((x, i) => a[i].mul(b[i])).reduce((m, n) => m.add(n));
  }

  static dotProduct2(a: number[], b: number[]): number {
    let sum = 0;
    const len = a.length;
    for (let i = 0; i < len; i++) {
      sum += a[i] * b[i];
    }
    return sum;
  }

  static dotProduct3(a: Decimal[], b: Decimal[]): Decimal {
    let result = new Decimal(0);
    const len = a.length;
    for (let i = 0; i < len; i++) {
      result = result.add(a[i].mul(b[i]));
    }
    return result;
  }

  static dotProduct4(a: Big[], b: Big[]): Big {
    let result = 0;
    const len = a.length;
    for (let i = 0; i < len; i++) {
      result = a[i].toNumber() * b[i].toNumber();
    }
    return new Big(result);
  }
}