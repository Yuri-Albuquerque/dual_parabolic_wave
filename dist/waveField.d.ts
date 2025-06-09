import { Point2D, WaveParams, SimulationConfig } from './types';
import { Parabola } from './parabola';
export declare class WaveField {
    private grid;
    private previousGrid;
    private velocityGrid;
    private boundaryMask;
    private config;
    private waveParams;
    private majorParabola;
    private minorParabola;
    private focusPoint;
    private time;
    constructor(config: SimulationConfig, waveParams: WaveParams, majorParabola: Parabola, minorParabola: Parabola, focusPoint: Point2D);
    private validateCFLCondition;
    private initializeGrids;
    private createBoundaryMask;
    private getGridCoordinates;
    private calculateSourceAmplitude;
    private applyWaveEquation;
    private applyBoundaryConditions;
    private addSourceExcitation;
    update(dt: number): void;
    getGrid(): number[][];
    getBoundaryMask(): boolean[][];
    getTime(): number;
    reset(): void;
}
//# sourceMappingURL=waveField.d.ts.map