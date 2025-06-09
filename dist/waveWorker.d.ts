interface WorkerMessage {
    type: 'init' | 'update' | 'reset' | 'setParams';
    data?: any;
}
interface WorkerResponse {
    type: 'initialized' | 'updated' | 'reset' | 'error';
    data?: any;
}
declare class WorkerWaveField {
    private grid;
    private previousGrid;
    private boundaryMask;
    private sourceGrid;
    private config;
    private waveParams;
    private focusPoint;
    private time;
    private focusI;
    private focusJ;
    constructor(config: any, waveParams: any, focusPoint: {
        x: number;
        y: number;
    });
    private initializeBoundaryMask;
    private calculateFocusPosition;
    private applyWaveEquation;
    private applyBoundaryConditions;
    private addSourceExcitation;
    update(dt: number): void;
    getGrid(): Float32Array;
    getBoundaryMask(): Uint8Array;
    getTime(): number;
    reset(): void;
}
declare let waveField: WorkerWaveField | null;
//# sourceMappingURL=waveWorker.d.ts.map