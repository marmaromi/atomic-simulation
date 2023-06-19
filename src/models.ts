export interface rVolumeModel {
  start: number;
  end: number;
  step: number;
}

export interface spheresResultsPropertiesModel {
  n: number;
  momX: string;
  momY: string;
  momZ: string;
  kinE: string;
  tCol: number;
  vSample: number;
  orderP: string;
  dv: number;
}

export interface spheresFinalResultsModel {
  'rVol': string;
  'sigma': string;
  'sum(time)': string;
  'sum(dv)': string;
  'pv0/kT': string;
}

