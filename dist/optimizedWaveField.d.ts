import { Point2D, WaveParams, SimulationConfig } from './types';
export declare class OptimizedWaveField {
    private worker;
    private grid;
    private boundaryMask;
    private config;
    private waveParams;
    private focusPoint;
    private time;
    private isInitialized;
    private updateCallbacks;
    constructor(config: SimulationConfig, waveParams: WaveParams, focusPoint: Point2D);
    private initializeWorker;
    private createWorkerBlob;
    onUpdate(callback: (grid: Float32Array, time: number) => void): void;
    update(dt: number, steps?: number): Promise<void>;
    reset(): Promise<void>;
    getGrid(): Float32Array;
    getBoundaryMask(): Uint8Array;
    getTime(): number;
    isReady(): boolean;
    destroy(): void;
}
//# sourceMappingURL=optimizedWaveField.d.ts.map