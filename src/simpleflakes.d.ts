declare module 'simpleflakes' {
  export class SimpleFlakeStruct {
    public timestamp: number;
    public randomBits: number;
    constructor(timestamp: number, randomBits: number);
  }

  export const simpleflake: (timestamp?: number, randomBits?: number, epoch?: number) => BigInt;
  export const SIMPLEFLAKE_EPOCH: number;
  export const binary: (val: string, padding: boolean) => string;
  export const parseSimpleflake: (val: string) => SimpleFlakeStruct;
}
