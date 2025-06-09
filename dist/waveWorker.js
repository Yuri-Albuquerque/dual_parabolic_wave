"use strict";
// Web Worker for wave field computation
// This runs in a separate thread to avoid blocking the main UI thread
// Simplified wave field implementation for worker
class WorkerWaveField {
    constructor(config, waveParams, focusPoint) {
        this.time = 0;
        this.focusI = 0;
        this.focusJ = 0;
        this.config = config;
        this.waveParams = waveParams;
        this.focusPoint = focusPoint;
        const gridSize = config.gridSize;
        this.grid = new Float32Array(gridSize * gridSize);
        this.previousGrid = new Float32Array(gridSize * gridSize);
        this.boundaryMask = new Uint8Array(gridSize * gridSize);
        this.sourceGrid = new Float32Array(gridSize * gridSize);
        this.initializeBoundaryMask();
        this.calculateFocusPosition();
    }
    initializeBoundaryMask() {
        const { gridSize, xRange, yRange } = this.config;
        const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
        const dy = (yRange[1] - yRange[0]) / (gridSize - 1);
        // Simplified parabola calculations for worker
        const majorDiameter = 20 * 25.4; // 508mm
        const majorFocus = 100;
        const minorDiameter = 100;
        const minorFocus = 50;
        for (let i = 0; i < gridSize; i++) {
            for (let j = 0; j < gridSize; j++) {
                const x = xRange[0] + j * dx;
                const y = yRange[1] - i * dy;
                // Major parabola (concave down): y = -x²/(4*f) + f
                const majorY = -(x * x) / (4 * majorFocus) + majorFocus;
                const insideMajor = y <= majorY && Math.abs(x) <= majorDiameter / 2;
                // Minor parabola (concave up): y = x²/(4*f) - f
                const minorY = (x * x) / (4 * minorFocus) - minorFocus;
                const outsideMinor = y >= minorY || Math.abs(x) > minorDiameter / 2;
                const index = i * gridSize + j;
                this.boundaryMask[index] = (insideMajor && outsideMinor) ? 1 : 0;
            }
        }
    }
    calculateFocusPosition() {
        const { gridSize, xRange, yRange } = this.config;
        const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
        const dy = (yRange[1] - yRange[0]) / (gridSize - 1);
        this.focusI = Math.round((yRange[1] - this.focusPoint.y) / dy);
        this.focusJ = Math.round((this.focusPoint.x - xRange[0]) / dx);
    }
    applyWaveEquation(dt) {
        const { gridSize, dampingFactor, xRange, yRange } = this.config;
        const c = this.waveParams.speed;
        const c2 = c * c;
        const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
        const dy = (yRange[1] - yRange[0]) / (gridSize - 1);
        const dx2 = dx * dx;
        const dy2 = dy * dy;
        // Wave equation coefficients based on C++ implementation
        const q0 = c * dt;
        const q1 = c2 * dt * dt;
        const q2 = (c * dt / dx) * (c * dt / dx);
        const q3 = (c * dt / dy) * (c * dt / dy);
        // Create new grid array
        const newGrid = new Float32Array(gridSize * gridSize);
        // Apply wave equation with optimized loop
        for (let i = 1; i < gridSize - 1; i++) {
            const iOffset = i * gridSize;
            const iPrevOffset = (i - 1) * gridSize;
            const iNextOffset = (i + 1) * gridSize;
            for (let j = 1; j < gridSize - 1; j++) {
                const index = iOffset + j;
                if (!this.boundaryMask[index]) {
                    newGrid[index] = 0;
                    continue;
                }
                // Finite difference operators
                const dxx = this.grid[iOffset + j - 1] - 2 * this.grid[index] + this.grid[iOffset + j + 1];
                const dyy = this.grid[iPrevOffset + j] - 2 * this.grid[index] + this.grid[iNextOffset + j];
                // Wave equation: u(t+dt) = (1 / (1 + damping*q0)) * 
                //                (-u(t-dt) + 2*u(t) + damping*q0*u(t-dt) + q1*source + q2*dxx + q3*dyy)
                const dampingTerm = dampingFactor * q0;
                const denominator = 1.0 + dampingTerm;
                newGrid[index] = (-this.previousGrid[index] +
                    2 * this.grid[index] +
                    dampingTerm * this.previousGrid[index] +
                    q1 * this.sourceGrid[index] +
                    q2 * dxx +
                    q3 * dyy) / denominator;
            }
        }
        // Update grids
        this.previousGrid.set(this.grid);
        this.grid.set(newGrid);
        // Apply boundary conditions
        this.applyBoundaryConditions();
    }
    applyBoundaryConditions() {
        const { gridSize, reflectionCoeff } = this.config;
        for (let i = 0; i < gridSize * gridSize; i++) {
            if (!this.boundaryMask[i]) {
                this.grid[i] = 0;
                this.previousGrid[i] = 0;
            }
        }
    }
    addSourceExcitation(time) {
        const { gridSize } = this.config;
        // Clear previous source
        this.sourceGrid.fill(0);
        if (this.focusI >= 0 && this.focusI < gridSize && this.focusJ >= 0 && this.focusJ < gridSize) {
            const index = this.focusI * gridSize + this.focusJ;
            if (this.boundaryMask[index]) {
                // Single pulse at initial time - Gaussian envelope with short duration
                const pulseWidth = 1.0 / this.waveParams.frequency; // One period width
                const pulseDuration = 2.0 * pulseWidth; // Duration of the pulse
                let amplitude = 0.0;
                if (time <= pulseDuration) {
                    // Gaussian envelope for smooth pulse
                    const gaussianWidth = pulseWidth / 3.0;
                    const envelope = Math.exp(-Math.pow(time - pulseWidth, 2) / (2 * gaussianWidth * gaussianWidth));
                    // Single frequency pulse
                    amplitude = this.waveParams.amplitude * 10.0 * envelope *
                        Math.sin(2 * Math.PI * this.waveParams.frequency * time);
                }
                this.sourceGrid[index] = amplitude;
                // Add to neighboring points for better excitation
                const neighbors = [
                    [-1, 0], [1, 0], [0, -1], [0, 1]
                ];
                for (const [di, dj] of neighbors) {
                    const ni = this.focusI + di;
                    const nj = this.focusJ + dj;
                    if (ni >= 0 && ni < gridSize && nj >= 0 && nj < gridSize) {
                        const neighborIndex = ni * gridSize + nj;
                        if (this.boundaryMask[neighborIndex]) {
                            this.sourceGrid[neighborIndex] = amplitude * 0.5;
                        }
                    }
                }
            }
        }
    }
    update(dt) {
        this.time += dt;
        this.addSourceExcitation(this.time);
        this.applyWaveEquation(dt);
    }
    getGrid() {
        return this.grid;
    }
    getBoundaryMask() {
        return this.boundaryMask;
    }
    getTime() {
        return this.time;
    }
    reset() {
        this.time = 0;
        this.grid.fill(0);
        this.previousGrid.fill(0);
        this.sourceGrid.fill(0);
    }
}
// Worker global variables
let waveField = null;
// Message handler
self.addEventListener('message', (event) => {
    const { type, data } = event.data;
    try {
        switch (type) {
            case 'init':
                const { config, waveParams, focusPoint } = data;
                waveField = new WorkerWaveField(config, waveParams, focusPoint);
                const response = {
                    type: 'initialized',
                    data: {
                        grid: waveField.getGrid(),
                        boundaryMask: waveField.getBoundaryMask(),
                        time: waveField.getTime()
                    }
                };
                self.postMessage(response);
                break;
            case 'update':
                if (waveField) {
                    const { dt, steps = 1 } = data;
                    // Perform multiple simulation steps
                    for (let i = 0; i < steps; i++) {
                        waveField.update(dt);
                    }
                    const updateResponse = {
                        type: 'updated',
                        data: {
                            grid: waveField.getGrid(),
                            time: waveField.getTime()
                        }
                    };
                    self.postMessage(updateResponse);
                }
                break;
            case 'reset':
                if (waveField) {
                    waveField.reset();
                    const resetResponse = {
                        type: 'reset',
                        data: {
                            grid: waveField.getGrid(),
                            time: waveField.getTime()
                        }
                    };
                    self.postMessage(resetResponse);
                }
                break;
            case 'setParams':
                if (waveField && data.waveParams) {
                    // Update wave parameters (requires recreation for simplicity)
                    const currentConfig = waveField['config'];
                    const currentFocus = waveField['focusPoint'];
                    waveField = new WorkerWaveField(currentConfig, data.waveParams, currentFocus);
                    const paramsResponse = {
                        type: 'initialized',
                        data: {
                            grid: waveField.getGrid(),
                            boundaryMask: waveField.getBoundaryMask(),
                            time: waveField.getTime()
                        }
                    };
                    self.postMessage(paramsResponse);
                }
                break;
        }
    }
    catch (error) {
        const errorResponse = {
            type: 'error',
            data: { message: error instanceof Error ? error.message : 'Unknown error' }
        };
        self.postMessage(errorResponse);
    }
});
//# sourceMappingURL=waveWorker.js.map