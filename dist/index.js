"use strict";
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __exportStar = (this && this.__exportStar) || function(m, exports) {
    for (var p in m) if (p !== "default" && !Object.prototype.hasOwnProperty.call(exports, p)) __createBinding(exports, m, p);
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.WaveField = exports.Parabola = exports.DualParabolicWaveSimulation = exports.ConsoleRenderer = void 0;
const simulation_1 = require("./simulation");
// For Node.js environment, we'll create a simple console-based visualization
class ConsoleRenderer {
    constructor() {
        this.intervalId = null;
        this.simulation = new simulation_1.DualParabolicWaveSimulation();
    }
    renderToConsole() {
        const grid = this.simulation['waveField'].getGrid();
        const boundaryMask = this.simulation['waveField'].getBoundaryMask();
        const gridSize = this.simulation['config'].gridSize;
        // Clear console
        console.clear();
        console.log('='.repeat(60));
        console.log('Dual Parabolic Wave Propagation Simulation');
        console.log(`Time: ${this.simulation['waveField'].getTime().toFixed(3)}s`);
        console.log(`Frequency: ${this.simulation.getWaveParams().frequency}Hz`);
        console.log('='.repeat(60));
        // Sample every few grid points for console display
        const step = Math.floor(gridSize / 40);
        let output = '';
        for (let i = 0; i < gridSize; i += step) {
            let line = '';
            for (let j = 0; j < gridSize; j += step) {
                const index = i * gridSize + j;
                if (boundaryMask[index]) {
                    const value = grid[index];
                    if (Math.abs(value) < 0.1) {
                        line += ' ';
                    }
                    else if (value > 0) {
                        line += '+';
                    }
                    else {
                        line += '-';
                    }
                }
                else {
                    line += '█'; // Boundary
                }
            }
            output += line + '\n';
        }
        console.log(output);
        // Print some statistics
        let maxAmplitude = 0;
        let activePoints = 0;
        for (let i = 0; i < gridSize; i++) {
            for (let j = 0; j < gridSize; j++) {
                const index = i * gridSize + j;
                if (boundaryMask[index]) {
                    activePoints++;
                    maxAmplitude = Math.max(maxAmplitude, Math.abs(grid[index]));
                }
            }
        }
        console.log(`Active grid points: ${activePoints}`);
        console.log(`Max amplitude: ${maxAmplitude.toFixed(4)}`);
        console.log('Legend: █=boundary, +=positive wave, -=negative wave, space=neutral');
    }
    startConsoleAnimation(intervalMs = 100) {
        this.intervalId = setInterval(() => {
            this.simulation['waveField'].update(0.001);
            this.renderToConsole();
        }, intervalMs);
    }
    stopConsoleAnimation() {
        if (this.intervalId) {
            clearInterval(this.intervalId);
            this.intervalId = null;
        }
    }
    getSimulation() {
        return this.simulation;
    }
}
exports.ConsoleRenderer = ConsoleRenderer;
// Export everything for use in browser or Node.js
var simulation_2 = require("./simulation");
Object.defineProperty(exports, "DualParabolicWaveSimulation", { enumerable: true, get: function () { return simulation_2.DualParabolicWaveSimulation; } });
var parabola_1 = require("./parabola");
Object.defineProperty(exports, "Parabola", { enumerable: true, get: function () { return parabola_1.Parabola; } });
var waveField_1 = require("./waveField");
Object.defineProperty(exports, "WaveField", { enumerable: true, get: function () { return waveField_1.WaveField; } });
__exportStar(require("./types"), exports);
// Main execution for Node.js
if (require.main === module) {
    console.log('Starting Dual Parabolic Wave Simulation...\n');
    console.log('Simulation Parameters:');
    console.log('- Major parabola: 20 inches (508mm) diameter, 100mm focus');
    console.log('- Minor parabola: 100mm diameter, 50mm focus');
    console.log('- Wave frequency: 1000Hz');
    console.log('- Sound speed: 343 m/s');
    console.log('\nPress Ctrl+C to stop the simulation\n');
    const renderer = new ConsoleRenderer();
    renderer.startConsoleAnimation(200);
    // Graceful shutdown
    process.on('SIGINT', () => {
        console.log('\nStopping simulation...');
        renderer.stopConsoleAnimation();
        process.exit(0);
    });
}
//# sourceMappingURL=index.js.map