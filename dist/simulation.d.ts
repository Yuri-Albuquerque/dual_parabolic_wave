import { Point2D, WaveParams } from './types';
import { Parabola } from './parabola';
export declare class DualParabolicWaveSimulation {
    private majorParabola;
    private minorParabola;
    private waveField;
    private focusPoint;
    private config;
    private waveParams;
    private animationId;
    private canvas;
    private ctx;
    private lastRenderTime;
    private renderInterval;
    private simulationStepsPerFrame;
    constructor();
    private setupParabolas;
    private setupWaveParameters;
    private setupSimulationConfig;
    private initializeWaveField;
    setupVisualization(canvasId: string): void;
    private drawParabolas;
    private drawWaveField;
    private render;
    private animate;
    start(): void;
    stop(): void;
    reset(): void;
    setFrequency(frequency: number): void;
    setAmplitude(amplitude: number): void;
    destroy(): void;
    getFocusPoint(): Point2D;
    getMajorParabola(): Parabola;
    getMinorParabola(): Parabola;
    getWaveParams(): WaveParams;
}
//# sourceMappingURL=simulation.d.ts.map