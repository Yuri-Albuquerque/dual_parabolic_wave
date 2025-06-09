import { DualParabolicWaveSimulation } from './simulation';
export declare class ConsoleRenderer {
    private simulation;
    private intervalId;
    constructor();
    renderToConsole(): void;
    startConsoleAnimation(intervalMs?: number): void;
    stopConsoleAnimation(): void;
    getSimulation(): DualParabolicWaveSimulation;
}
export { DualParabolicWaveSimulation } from './simulation';
export { Parabola } from './parabola';
export { WaveField } from './waveField';
export * from './types';
//# sourceMappingURL=index.d.ts.map